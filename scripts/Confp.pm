#!/usr/bin/env perl

package Confp;

use strict;
use warnings;

sub new {
	my $proto = shift;
	my $class = ref($proto) || $proto;
	my $self = {};
	
	$self->{PATH} = "";
	$self->{CONFIG} = {};
	bless($self, $class);

	if (@_) {
		if (!$self->read(shift)) {
			return undef;
		}
	}

	return $self;
}

sub read {
	my $self = shift;
	my $path = shift;

	open(CONFFILE, "<", $path) || return undef;

	$self->{PATH} = $path;

	my $n = 1;
	while (<CONFFILE>) {
		if (!$self->parse($_)) {
			print "Syntax error at line ", $n;
			return 0;
		}
		$n++;
	}
	
	close(CONFFILE);
	
	return 1;
}

sub trim {
	my $string = shift;
	$string =~ s/^\s+//;
	$string =~ s/\s+$//;
	return $string;
}

sub ltrim {
	my $string = shift;
	$string =~ s/^\s+//;
	return $string;
}

sub rtrim {
	my $string = shift;
	$string =~ s/\s+$//;
	return $string;
}

sub parse { 
	my $self = shift;
	my $line = trim(shift);

	if (substr ($line, 0, 1) =~ /#/) {
		return 1;
	} elsif ($line =~ /\w+\s*=\s*(\".*\"|\'.*\'|[^\s\"\']+)/) {

		my ($key, $value) = split(/=/, $line, 2);

		$key = trim($key);
		$value = trim($value);
		
		my $rsqi = rindex($value, "'");
		my $lsqi = index($value, "'");
		my $rdqi = rindex($value, "\"");
		my $ldqi = index($value, "\"");
		my $ci = index($value, "#");

		if ($lsqi == 0) {
			if ($rsqi < 0) {
				return 0;
			}
				
			$value = substr($value, $lsqi + 1, $rsqi - $lsqi - 1);
		} elsif ($ldqi == 0) {
			if ($ldqi < 0) {
				return 0;
			}

			$value = substr($value, $ldqi + 1, $rdqi - $ldqi - 1);
		} elsif ($ci >= 0) {
			if ($ci == 0) {
				return 0;
			}

			$value = substr($value, 0, $ci - 1);
		}

		$self->{CONFIG}{$key} = $value;

		return 1;
	} elsif (!$line) {
		return 1;
	} else {
		return 0;
	}
}


sub get { 
	my $self = shift;
	my $key  = shift;

	return $self->{CONFIG}{$key};
}

1;
