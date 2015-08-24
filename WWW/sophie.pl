#!/usr/bin/perl
#*************************************************************************
#
#   Program:    sophie
#   File:       sophie.pl
#   
#   Version:    V1.0
#   Date:       12.06.97
#   Function:   HTTP query server for finding human subgroups
#   
#   Copyright:  (c) UCL / Dr. Andrew C. R. Martin 1997
#   Author:     Dr. Andrew C. R. Martin
#   Address:    Biomolecular Structure & Modelling Unit,
#               Department of Biochemistry & Molecular Biology,
#               University College,
#               Gower Street,
#               London.
#               WC1E 6BT.
#   EMail:      andrew@bioinf.org.uk
#               
#*************************************************************************
#
#   This program is not in the public domain, but it may be copied
#   according to the conditions laid out in the accompanying file
#   COPYING.DOC
#
#   The code may be modified as required, but any modifications must be
#   documented so that the person responsible can be identified. If 
#   someone else breaks this code, I don't want to be blamed for code 
#   that does not work! 
#
#   The code may not be sold commercially or included as part of a 
#   commercial product except as described in the file COPYING.DOC.
#
#*************************************************************************
#
#   Description:
#   ============
#
#*************************************************************************
#
#   Usage:
#   ======
#
#*************************************************************************
#
#   Revision History:
#   =================
#   V1.0  12.06.97 Original    By: ACRM
#
#*************************************************************************
# Main Program
# ------------
# Sets names of external programs and environment variables they require.
# Gets and parses the HTTP form input information. Build a sequence file
# and tests it using the sophie program.
# 
# 12.06.97 Original    By: ACRM
# 
# The user to receive messages when a test has been run
$MailRecipient = "amartin";
# Don't send mail if query came from this host
$NoMailHost    = "bruno";

# External Programs
$bindir    = "/user/amartin/sophie";
$sophie    = "$bindir/sophie";

# Get environment variables from the HTTP form
$method     = $ENV{'REQUEST_METHOD'};
$type       = $ENV{'CONTENT_TYPE'};
$length     = $ENV{'CONTENT_LENGTH'};
$RemoteHost = $ENV{'REMOTE_HOST'};

# Turn on flushing
$| = 1;

# Build the filename for the temporary files
$seqfile = "/tmp/sophie.in.$$";

# Output MIME header information
print "Content-type: text/html\n\n";

# Check the method and type are correct
if($method ne "POST")
{
    print "<p><h2>Error:</h2>\n";
    print "This script may only be referenced using the http POST method.\n";
    print 'See <a href="http://www.ncsa.uiuc.edu/SDG/Software/Mosaic/Docs/fill-out-forms/overview.html">forms overview</a> for details';
    print "</p>\n\n";
    print "<p>Supplied method was: $method</p>\n\n";
    exit 1;
}
if($type ne "application/x-www-form-urlencoded")
{
    print "<p><h2>Error:</h2>\n";
    print "This script can only be used to decode http FORM results.</p>\n";
    exit 1;
}

# Parse the input
%CgiData = &ParseInput;
$address = $CgiData{'address'};

# Check that the sequences have both been specified
if(! &CheckInput(%CgiData))
{
    print "<p><h2>Error:</h2>\n";
    print "You must provide a sequence!</p>\n";
    exit 1;
}

# Create a sequence file
&BuildSequenceFile($seqfile, %CgiData);

# Print a header
print "<h1>Human Subgroup Assignment Results</h1>\n";
print "<p>If no results appear below, please send me EMail with your\n";
print "sequence as an error has occurred!</p>\n\n";

# Create the Kabat numbered sequence file
system("$sophie $seqfile");



# Remove the temporary files
unlink($seqfile);

# Send an EMail message to say that a sequence analysis has been done
if($RemoteHost ne $NoMailHost)
{
    &SendMail($MailRecipient, $address, $RemoteHost,
              $CgiData{'sequence'});
}

## END OF MAIN PROGRAM


#*************************************************************************
# %cgidata ParseInput()
# ---------------------
# Routine to parse the input and build an associative array of data
# items and their contents.
# 00.00.00 Original   By: Alex Michie
# 14.03.95 Tidied up, added comments and changed to use split() By: ACRM
# 20.03.95 Modified so = signs will work properly in the data   By: ACRM
#
sub ParseInput 
{
    local($req,@data,$i,$key,$val,%output);

    # Place the input data into the $req string
    if ($ENV{'REQUEST_METHOD'} eq "GET") 
    {
        $req = $ENV{'QUERY_STRING'};
    }
    elsif ($ENV{'REQUEST_METHOD'} eq "POST") 
    {
        read(STDIN,$req,$ENV{'CONTENT_LENGTH'});  
    }
    
    # Replace + signs with spaces
    $req =~ s/\++/ /g;  
           
    # Split at & signs into the data array
    @data = split(/&/,$req);     

    # For each item in the data array
    foreach $i (0..$#data)
    {
        # Strip leading and trailing spaces
        $data[$i] =~ s/ +$//;    
        $data[$i] =~ s/= +/=/;                   

        # Split into key and value
        ($key, $val) = split(/=/,$data[$i]);
        # Filter hex to something sensible
        $key =~ s/%(..)/pack("c",hex($1))/ge;  
        $val =~ s/%(..)/pack("c",hex($1))/ge;  
        # Dump into an associative array for return
        $output{$key} = $val;
    }

    # Return the associative array
    return(%output);
}


#*************************************************************************
# void BuildSequenceFile($seqfile, %CgiData)
# ------------------------------------------
# Build a sequence file for submitting to the test program
# 12.06.97 Original    By: ACRM
# 
sub BuildSequenceFile
{
    local($seqfile, %CgiData) = @_;

    # Open a file to store the sequence
    open(SEQFILE,">$seqfile") || die "Unable to write sequence file\n";

    # Place header into file
    print SEQFILE ">P1;ABSEQ\n";
    print SEQFILE "Antibody sequence - $CgiData{'address'}\n";

    # Print the light chain and a * if necessary
    print SEQFILE "$CgiData{'sequence'}";
    if($CgiData{'sequence'} =~ /\*/)
    {
        print SEQFILE "\n";
    }
    else
    {
        print SEQFILE "*\n";
    }
    
    # Close the file
    close(SEQFILE);
}


#*************************************************************************
# BOOL CheckInput
# ---------------
# Checks that a sequence has been supplied
# 12.06.97 Original    By: ACRM
# 
sub CheckInput
{
    local(%data) = @_;

    if($data{'sequence'} ne "")
    {
        return 1;
    }

    return 0;
}


#*************************************************************************
# void SendMail($recipient, $address, $RemoteHost, $light, $heavy)
# ----------------------------------------------------------------
# Sends a mail message to $recipient to say that $address has done a
# search.
# 14.03.95 Original    By: ACRM
# 15.03.95 Modified to display sequences
# 
sub SendMail
{
    local($recipient, $address, $RemoteHost, $sequence) = @_;
    local($message);

    $message = <<"EOF";
Subject: Subgroup assignment
To: $recipient

A Subgroup assignment has been run (from node $RemoteHost) by $address 

Sequence was: 
$sequence

EOF

    open(MAIL, "|/usr/lib/sendmail -t") || die " ";
    print MAIL $message;
    close(MAIL);
}


