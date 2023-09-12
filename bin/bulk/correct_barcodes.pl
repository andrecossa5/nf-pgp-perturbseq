use strict;

my $whitelist_line;
my @whitelist_line_split;
my %hash_table;
my $GBC;

## read whitelist
open(IN1, "<", $ARGV[0]) or die $!;
while($whitelist_line = <IN1>)
{
	## remove trailing newline character
	chomp($whitelist_line);
	## split line by space
	@whitelist_line_split = split(/\s+/, $whitelist_line);
	## add first column (degenerated GBC) to hash table
	$hash_table{$whitelist_line_split[0]} = $whitelist_line_split[1];
}

## read GBCs
open(IN2, "gunzip -c $ARGV[1] |") or die $!;
while($GBC = <IN2>)
{
	## remove trailing newline character
	chomp($GBC);
	## if GBC is present in hash table, retrieve correct GBC, else print the
	## non-corrected GBC
	if (exists $hash_table{$GBC})
	{
		print("$hash_table{$GBC}\n");
	}
	else
	{
		print("$GBC\n");
	}
}
