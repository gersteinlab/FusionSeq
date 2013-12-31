#!/usr/bin/env perl

use strict;
use warnings;

use FindBin;

use lib "$FindBin::Bin";

use Confp;
use Cwd;

use constant CGIS => "geneFusions_cgi showDetails_cgi seqViz_cgi findFusionPartner_cgi";

my $path = getcwd() . "/default.fusionseqrc";
my $config = Confp->new($path);

die("Cannot open default.fusionseqrc") if !$config;

my $cmd = "mv " . CGIS . " " . $config->get('WEB_DATA_DIR');

#print $cmd, "\n";
system($cmd);
