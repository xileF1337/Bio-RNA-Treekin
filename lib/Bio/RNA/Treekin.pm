package Bio::RNA::Treekin;

use 5.006;
use strict;
use warnings;

our $VERSION = '0.01';


# Stores a data from a single row of the treekin file, i.e. the populations of
# all minima at a given time point.
package Bio::RNA::Treekin::PopulationDataRecord {
    use Moose;
    use MooseX::StrictConstructor;
    use namespace::autoclean;

    use autodie qw(:all);
    use Scalar::Util qw(reftype);
    use   List::Util qw(    all);

    use overload '""' => \&stringify;

    has 'time'  => (is => 'ro', required => 1);

    has '_populations' => (
        is       => 'ro',
        required => 1,
        init_arg => 'populations',
    );

    sub clone {
        my $self = shift;
        my $clone = __PACKAGE__->new(
            time        => $self->time,
            populations => [ $self->populations ],
        );
        return $clone;
    }

    sub min_count {
        my $self = shift;
        my $min_count = @{ $self->_populations };   # number of data points

        return $min_count;
    }

    # Use to adjust the min count, e.g. when the passed data array was
    # constructed before the number of minima was known. It may not be
    # shrinked as data might be lost.
    sub set_min_count {
        my ($self, $new_min_count) = @_;

        my $current_min_count = @{ $self->_populations };
        confess 'Can only increase min_count'
            if  $current_min_count > $new_min_count;

        # Set additional states to population of 0.
        for my $i ( $current_min_count..($new_min_count-1) ) {
            $self->_populations->[$i] = 0.;
        }

        return;
    }

    # Return populations of all mins. Use of_min() instead to get the
    # population of a specific min.
    # Returns a list of all minima's populations.
    sub populations {
        my $self = shift;
        return @{ $self->_populations };
    }

    # Get population for the given minimum.
    sub of_min {
        my ($self, $min) = @_;

        confess "Minimum $min is out of bounds"
            if $min < 1 or $min > $self->min_count;

        # Minimum 1 is the first one (index 0)
        my $population = $self->_populations->[$min-1];
        return $population;
    }

    # Transform (reorder and resize) the population data according to a given
    # mapping and resize to a given minimum count. All minima that are not
    # mapped to a new position are discarded (replaced by zero).
    # NOTE: Ensure that no two minima are mapped to the same position or crap
    # will happen.
    # Arguments:
    #   maps_to_min_ref: Hash ref that specifies for each kept minimum (key) to
    #       which new minimum (value) it is supposed to be mapped.
    #   new_min_count: New size (number of mins) of the record after the
    #       transformation.
    # Void.
    sub transform {
        my ($self, $maps_to_min_ref, $new_min_count) = @_;

        my @new_pops = (0) x $new_min_count;    # new array with right size
        # my @target_indices
        #     = map { $maps_to_min_ref->{$_} - 1 } 1..$self->min_count;
        my @source_mins    = grep { defined $maps_to_min_ref->{$_} }
                                  1..$self->min_count;      # filter unmapped
        my @source_indices = map  { $_ - 1 } @source_mins;
        my @target_indices = map  { $maps_to_min_ref->{$_} - 1 } @source_mins;

        # Sanity check.
        confess "Cannot reorder as some minima are not mapped correctly"
            unless all {$_ >= 0 and $_ < $new_min_count} @target_indices;

        # Copy population data to the correct positions and overwrite original.
        @new_pops[@target_indices] = @{$self->_populations}[@source_indices];
        # @new_pops[@target_indices] = @{ $self->_populations };
        @{ $self->_populations }   = @new_pops;

        return;
    }

    sub _parse_population_data_line {
        my ($self, $population_data_line) = @_;

        my ($time, @populations) = split /\s+/, $population_data_line;

        my @args = (
            time        => $time,
            populations => \@populations,
        );

        return @args;
    }

    around BUILDARGS => sub {
        my $orig = shift;
        my $class = shift;

        return $class->$orig(@_) if @_ != 1 or reftype $_[0];

        # We have a population data line here.
        my $population_data_line = shift;
        my @args
            = $class->_parse_population_data_line($population_data_line);

        return $class->$orig(@args);
    };

    sub stringify {
        my $self = shift;

        # Format data like treekin C code
        my $self_as_string = sprintf "%22.20e ", $self->time;
        $self_as_string
            .= join q{ }, map {sprintf "%e", $_} @{ $self->_populations };

        # There seems to be a trailing space in the treekin C code
        # (printf "%e ") but there is none in the treekin simulator output.
        # $self_as_string .= q{ };

        return $self_as_string;
    }

    __PACKAGE__->meta->make_immutable;
} # End of Bio::RNA::Treekin::PopulationDataRecord

package Bio::RNA::Treekin::Record {
    use Moose;
    use MooseX::StrictConstructor;
    use namespace::autoclean;

    use autodie qw(:all);
    use Scalar::Util qw(reftype openhandle);
    use List::Util qw(first pairmap max uniqnum all);
    use Carp qw(croak);

    use overload '""' => \&stringify;

    has '_population_data'  => (
        is       => 'ro',
        required => 1,
        init_arg => 'population_data',
    );

    has 'date'              => (is => 'ro', required => 1);
    has 'sequence'          => (is => 'ro', required => 1);
    has 'method'            => (is => 'ro', required => 1);
    has 'start_time'        => (is => 'ro', required => 1);
    has 'stop_time'         => (is => 'ro', required => 1);
    has 'temperature'       => (is => 'ro', required => 1);
    has 'basename'          => (is => 'ro', required => 1);
    has 'time_increment'    => (is => 'ro', required => 1);
    has 'degeneracy'        => (is => 'ro', required => 1);
    has 'absorbing_state'   => (is => 'ro', required => 1);
    has 'states_limit'      => (is => 'ro', required => 1);

    # Add optional attributes including predicate.
    has $_ => (
                   is        => 'ro',
                   required  => 0,
                   predicate => "has_$_",
              )
        foreach qw(
                     info
                     init_population
                     rates_file
                     file_index
                     cmd
                );

    # Get number of population data rows stored.
    sub population_data_count {
        my ($self) = @_;

        my $data_count = @{ $self->_population_data };
        return $data_count;
    }

    # Number of states / minima in this simulation.
    # Get number of mins in the first population record; it should be the
    # same for all records.
    sub min_count {
        my $self = shift;

        my $first_pop = $self->population(0);
        confess 'min_count: no population data present'
            unless defined $first_pop;

        my $min_count = $first_pop->min_count;

        return $min_count;
    }

    # Keep only the population data for the selected minima, remove all other.
    # Will NOT rescale populations, so they may no longer sum up to 1.
    # Arguments:
    #   mins: List of mins to keep. Will be sorted and uniq'ed (cf. splice()).
    # Returns the return value of splice().
    sub keep_mins {
        my ($self, @kept_mins) = @_;
        @kept_mins = uniqnum sort {$a <=> $b} @kept_mins;   # sort / uniq'ify
        return $self->splice_mins(@kept_mins);
    }

    # Keep only the population data for the selected minima, remove all other.
    # May duplicate and re-order.
    #   mins: List of mins to keep. Will be used as is.
    # Returns itself.
    sub splice_mins {
        my ($self, @kept_mins) = @_;

        my $min_count = $self->min_count;
        confess 'Cannot splice, minimum out of bounds'
            unless all {$_ >= 1 and $_ <= $min_count} @kept_mins;

        # Directly update raw population data here instead of doing tons of
        # calls passing the same min array.
        my @kept_indices = map {$_ - 1} @kept_mins;
        for my $pop_data (@{$self->_population_data}) { # each point in time
            my $raw_pop_data = $pop_data->_populations;
            @{$raw_pop_data} = @{$raw_pop_data}[@kept_indices];
        }

        return $self;
    }

    # Get the maximal population for the given minimum over all time points.
    sub max_pop_of_min {
        my ($self, $min) = @_;
        my $max_pop = '-Inf';
        for my $pop_data (@{$self->_population_data}) { # each point in time
            $max_pop = max $max_pop, $pop_data->of_min($min);   # update max
        }
        return $max_pop;
    }

    # For a given minimum, return all population values in chronological
    # order.
    # Arguments:
    #   min: Minimum for which to collect the population data.
    # Returns a list of population values in chronological order.
    sub pops_of_min {
        my ($self, $min) = @_;

        my @pops_of_min = map { $_->of_min($min) } $self->populations;

        return @pops_of_min;
    }

    # Final population data record, i.e. the result of the simulation.
    sub final_population {
        my ($self) = @_;

        my $final_population_data
            = $self->population($self->population_data_count - 1);

        return $final_population_data;
    }

    # Get the i-th population data record (0-based indexing).
    sub population {
        my ($self, $i) = @_;

        my $population_record = $self->_population_data->[$i];
        return $population_record;
    }

    # Return all population data records.
    sub populations {
        return @{ $_[0]->_population_data };
    }

    # Add a new minimum with all-zero entries. Data can then be appended to
    # this new min.
    # Returns the index of the new minimum.
    sub add_min {
        my $self = shift;
        my $new_min_count = $self->min_count + 1;

        # Increase the min count of all population data records by one.
        $_->set_min_count($new_min_count)
            foreach @{ $self->_population_data }, $self->init_population;

        return $new_min_count;              # count == highest index
    }

    # Given a key--value list of minima (key) and population data list refs
    # (value), append the data points to the respective minimum, filling other
    # columns with zeroes.
    # The passed population data objects are modified.
    sub append_pop_data {
        my ($self, $pop_data_ref, $merged_to_min_ref) = @_;

        my $min_count = $self->min_count;
        $_->transform($merged_to_min_ref, $min_count) foreach @$pop_data_ref;

        push @{ $self->_population_data }, @$pop_data_ref;

        return;
    }

    # Decode a single header line into a key and a value, which are returned.
    sub _get_header_line_key_value {
        my ($class, $header_line) = @_;

        # key and value separated by first ':' (match non-greedy!)
        my ($key, $value) = $header_line =~ m{ ^ ( .+? ) : [ ] ( .* ) $ }x;

        confess "Invalid key in header line:\n$header_line"
            unless defined $key;

        # Convert key to lower case and replace spaces by underscores.
        $key = (lc $key) =~ s/\s+/_/gr;

        return ($key, $value);
    }

    # Decode the initial population from the Treekin command line. The
    # population is given as multiple --p0 a=x switches, where a is the state
    # index and x is the fraction of population initially present in this
    # state.
    # Arguments:
    #   command: the command line string used to call treekin
    # Returns a hash ref containing the initial population of each state a at
    # position a (1-based).
    sub _parse_init_population_from_cmd {
        my ($class, $command) = @_;

        my @command_parts = split /\s+/, $command;

        # Extract the initial population strings given as (multiple) arguments
        # --p0 to Treekin from the Treekin command.
        my @init_population_strings;
        while (@command_parts) {
            if (shift @command_parts eq '--p0') {
                # Next value should be a population value.
                confess 'No argument following a --p0 switch'
                    unless @command_parts;
                push @init_population_strings, shift @command_parts;
            }
        }

        # Store population of state i in index i-1.
        my @init_population;
        foreach my $init_population_string (@init_population_strings) {
            my ($state, $population) = split /=/, $init_population_string;
            $init_population[$state-1] = $population;
        }

        # If no population was specified on the cmd line, init 100% in state 1
        $init_population[0] = 1 unless @init_population_strings;

        # Set undefined states to zero.
        $_ //= 0. foreach @init_population;

        my $init_population_record
            = Bio::RNA::Treekin::PopulationDataRecord->new(
                time        => 0,
                populations => \@init_population,
        );
        return $init_population_record;
    }

    sub _parse_header_lines {
        my ($class, $header_lines_ref) = @_;

        my @header_args;
        foreach my $line (@$header_lines_ref) {
            my ($key, $value) = $class->_get_header_line_key_value($line);

            # Implement special handling for certain keys.
            if ($key eq 'rates_file') {
                # remove (#index) from file name and store the value
                my ($file_name, $file_index)
                    = $value =~ m{ ^ (.+) [ ] [(] [#] (\d+) [)] $ }x;
                push @header_args, (
                    rates_file  => $file_name,
                    file_index  => $file_index,
                );
            }
            elsif ($key eq 'cmd') {
                # Extract initial population from Treekin command.
                my $init_population_ref
                    = $class->_parse_init_population_from_cmd($value);
                push @header_args, (
                    cmd             => $value,
                    init_population => $init_population_ref,
                );
            }
            else {
                # For the rest, just push key and value as constructor args.
                push @header_args, ($key => $value);
            }
        }
        return @header_args;
    }


    # Read all lines from the given handle and separate it into header lines
    # and data lines.
    sub _read_record_lines {
        my ($class, $record_handle) = @_;

        # Separate lines into header and population data.  All header lines
        # begin with a '# ' (remove it!)
        my ($current_line, @header_lines);
        while (defined ($current_line = <$record_handle>)) {
            next if $current_line =~ /^@/;  # drop xmgrace annotations
            # header lines start with '# ', remove it
            last unless $current_line =~ s/^# //;
            push @header_lines, $current_line;
        }
        chomp @header_lines;

        #say STDERR '<', @header_lines, '>' and
        confess 'Unexpected end of record while parsing header'
            unless defined $current_line;

        my @population_data_lines = ($current_line);
        while (defined ($current_line = <$record_handle>)) {
            push @population_data_lines, $current_line;
        }
        chomp @population_data_lines;

        return \@header_lines, \@population_data_lines;
    }

    sub _parse_population_data_lines {
        my ($class, $population_data_lines_ref) = @_;

        my @population_data
            = map { Bio::RNA::Treekin::PopulationDataRecord->new($_) }
                  @$population_data_lines_ref
                  ;

        return (population_data => \@population_data);
    }

    around BUILDARGS => sub {
        my $orig  = shift;
        my $class = shift;

        # Call original constructor if passed more than one arg.
        return $class->$orig(@_) unless @_ == 1;

        # Retrive file handle or pass on hash ref to constructor.
        my $record_handle;
        if (reftype $_[0]) {
            if (reftype $_[0] eq reftype {}) {          # arg hash passed,
                return $class->$orig(@_);               # pass on as is
            }
            elsif (reftype $_[0] eq reftype \*STDIN) {  # file handle passed
                $record_handle = shift;
            }
            else {
                croak 'Invalid ref type passed to constructor';
            }
        }
        else {                                          # file name passed
            my $record_file = shift;
            open $record_handle, '<', $record_file;
        }

        # Read in file.
        my ($header_lines_ref, $population_data_lines_ref)
            = $class->_read_record_lines($record_handle);

        # Parse file.
        my @header_args = $class->_parse_header_lines($header_lines_ref);
        my @data_args
            = $class->_parse_population_data_lines($population_data_lines_ref);

        my %args = (@header_args, @data_args);
        return $class->$orig(\%args);
    };

    sub BUILD {
        my $self = shift;

        # Force construction despite lazyness.
        $self->min_count;

        # Adjust min count of initial population as it was not known when
        # initial values were extracted from treekin cmd.
        $self->init_population->set_min_count( $self->min_count )
            if $self->has_init_population;
    }

    sub stringify {
        my $self = shift;

        # Format header line value of rates file entry.
        my $make_rates_file_val = sub {
            $self->rates_file . ' (#' . $self->file_index . ')';
        };

        # Header
        my @header_entries = (
            $self->has_rates_file ? ('Rates file' => $make_rates_file_val->()) : (),
            $self->has_info       ? ('Info'       => $self->info)   : (),
            $self->has_cmd        ? ('Cmd'        => $self->cmd)    : (),
            'Date'            => $self->date,
            'Sequence'        => $self->sequence,
            'Method'          => $self->method,
            'Start time'      => $self->start_time,
            'Stop time'       => $self->stop_time,
            'Temperature'     => $self->temperature,
            'Basename'        => $self->basename,
            'Time increment'  => $self->time_increment,
            'Degeneracy'      => $self->degeneracy,
            'Absorbing state' => $self->absorbing_state,
            'States limit'    => $self->states_limit,
        );

        my $header_str = join "\n", pairmap { "# $a: $b" } @header_entries;

        # Population data
        my $population_str
            = join "\n", map { "$_" } @{ $self->_population_data };

        my $self_as_str = $header_str . "\n" . $population_str;
        return $self_as_str;
    }

    __PACKAGE__->meta->make_immutable;
};  # End of Bio::RNA::Treekin::Record


package IO::File::RecordStream {
    use Moose;
    use MooseX::StrictConstructor;
    use namespace::autoclean;

    use autodie qw(:all);
    use Scalar::Util qw(reftype openhandle);
    #use List::Util qw(any);

    use IO::Lines;

    has 'file_name' => (
        is        => 'ro',
        isa       => 'Str',
        predicate => 'has_file_name',
    );

    has 'file_handle' => (
        is          => 'ro',
        isa         => 'FileHandle',
        builder     => '_build_file_handle',
        lazy        => 1,
    );

    has 'end_reached' => (
        is          => 'ro',
        default     => 0,
        init_arg    => undef,
        writer      => '_end_reached',           # private writer
    );

    # A regexp matching the separator line used to separate individual records
    has 'match_separator' => (
        is          => 'ro',
        isa         => 'RegexpRef',
        required    => 1,
    );

    # A code ref that can be passed a ref to the array containing the read
    # lines and that makes a new record object from it.
    has '_record_factory' => (              # keep the ref private
        is          => 'ro',
        isa         => 'CodeRef',
        init_arg    => 'record_factory',
        required    => 1,
    );

    # Allow various calling styles of the constructor:
    # new(file_handle): pass file handle to read data from
    # new(file_name):   pass file name of file to read data from
    around BUILDARGS => sub {
        my $orig  = shift;
        my $class = shift;

        return $class->$orig(@_) unless @_ == 1;  # no special handling

        # Check if we got a file name or handle for multi-record input file.
        my @constructor_args;
        if (not reftype $_[0]) {                    # file name given
            my $input_file_name = shift;
            push @constructor_args, (file_name => $input_file_name);
        }
        elsif (reftype $_[0] eq reftype \*STDIN) {  # file handle given
            my $input_file_handle = shift;
            push @constructor_args, (file_handle => $input_file_handle);
        }
        else {                                      # no file name / handle
            return $class->$orig(@_);
        }

        return $class->$orig(@constructor_args);
    };


    sub BUILD {
        my $self = shift;

        confess 'The value of file_handle does not seem to be an open handle'
            unless openhandle $self->file_handle;

        # If the input file is empty, set end_reached immediately.
        $self->_end_reached(1) if eof $self->file_handle;

        return;
    }

    # Open file handle from file name if no handle was passed. Die if we cant.
    sub _build_file_handle {
        my $self = shift;

        confess 'Cannot build file handle unless a file name was specified'
            unless $self->has_file_name;

        open my $input_file_handle, '<', $self->file_name;
        return $input_file_handle;
    }

    # Returns chomped lines.
    sub _read_next_record {
        my $self = shift;

        my $record_file_handle = $self->file_handle;
        my @record_lines;
        while (<$record_file_handle>) {
            my $line = $_;
            chomp $line;

            if ($line =~ $self->match_separator) {            # end of record
                # Drop separator line and return current collection of lines.
                # Also test if file ends in a separator.
                $self->_end_reached(1) if eof $record_file_handle;
                return \@record_lines
            }
            else {
                # Store lines until end of record
                push @record_lines, $line;
            }
        }

        $self->_end_reached(1);                     # file has been read
        return \@record_lines;
    }

    # Get the next record from the multi-record file.
    sub next {
        my $self = shift;

        # Are there any more entries?
        return if $self->end_reached;

        # Read lines of next record.
        my $record_lines_ref = $self->_read_next_record;

        # Construct record object using factory
        my $record_array_handle = IO::Lines->new($record_lines_ref);
        my $record = $self->_record_factory->($record_array_handle);

        return $record;
    }

    __PACKAGE__->meta->make_immutable;
}; # End of IO::File::RecordStream


package Bio::RNA::Treekin::MultiRecord {
    use Moose;
    use MooseX::StrictConstructor;
    use namespace::autoclean;

    use autodie qw(:all);

    extends 'IO::File::RecordStream';

    has '+_record_factory' => (     # + means overwrite inherited attribute
        is       => 'ro',
        init_arg => undef,
        default  => sub { return sub { Bio::RNA::Treekin::Record->new(@_); } },
    );

    has '+match_separator' => (     # + means overwrite inherited attribute
        is       => 'ro',
        init_arg => undef,
        default  => sub {
                         qr{ ^ & $ }x  # match a line consisting of single '&'
                    },
    );

    __PACKAGE__->meta->make_immutable;
} # End of Bio::RNA::Treekin::MultiRecord


=head1 NAME

Bio::RNA::Treekin - Classes for working with I<Treekin> output.

=head1 VERSION

Version 0.01

=cut


=head1 SYNOPSIS

This module provides auxiliary classes to parse, query and print the data
generated by the reaction kinetics simulation tool I<Treekin>, as well as
multi-simulations as generated by the I<BarMap> simulator.

Perhaps a little code snippet.

    use Bio::RNA::Treekin;

    # Read multi-record file (records separated by single '&')
    my $treekin_records = Bio::RNA::Treekin::MultiRecord->new($treekin_file);

    # Iterate over single records.
    while (defined (my $next_treekin_record = $treekin_records->next)) {
        # Print current file name.
        print 'Treekin simulation for ', $next_treekin_record->rates_file, " \n";

        # Print population of minimum 42 at the end of the simulation
        my  $init_pop = $current_treekin_record->init_population;
        my $final_pop = $current_treekin_record->final_population;
        print 'min 42: ', $init_pop->of_min(42), ' => ',
              $final_pop->of_min(42), "\n";
    }


=head1 METHODS OF C<Bio::RNA::Treekin::MultiRecord>

Reads multiple I<Treekin> records from a single file. The individual records
are separated by a line containing a single '&'.

=head2 C<new()>

Constructor. Pass a file name or a file handle to the I<Treekin> multi-record
file.

=head2 C<next()>

Returns the next I<Treekin> record, i.e. an object of type
C<Bio::RNA::Treekin::Record>.


=head1 METHODS OF C<Bio::RNA::Treekin::Record>


=head2 C<init_population()>

Get initial population data record, i.e. an object of type
C<Bio::RNA::PopulationDataRecord>. This is equal to the values given to
I<Treekin> via the C<--p0> command line arguments.

=head2 C<population_data_count()>

Number of population data records, i.e. number of time points for which
population data has been computed.

=head2 C<min_count()>

Number of states / minima in this simulation.

=head2 C<keep_mins(@mins)>

Keep only the population data for the selected minima, remove all other.
Will NOT rescale populations, so they may no longer sum up to 1.

Arguments:

=over

=item @mins

List of mins to keep. Will be sorted and uniq'ed (cf. C<splice()>).

=back

Returns the return value of C<splice()>.

=head2 C<splice_mins()>

Keep only the population data for the selected minima, remove all other.
May duplicate and re-order.

Arguments:

=over

=item @mins

List of mins to keep. Will be used as is.

=back

Returns itself.

=head2 C<max_pop_of_min()>

Get the maximal population for the given minimum over all time points.

=head2 C<final_population()>

Final population data record, i.e. the result of the simulation.

=head2 C<population()>

Get the i-th population data record (0-based indexing).

=head2 C<stringify()>

Converts the object back into a string. Also called by stringification
overloading in double quotes.

=head2 Querying header information

There are accessor methods for following attributes given in the header.

=over

=item *
rates_file

=item *

file_index

=item *

cmd

=item *

date

=item *

sequence

=item *

method

=item *

start_time

=item *

stop_time

=item *

temperature

=item *

basename

=item *

time_increment

=item *

degeneracy

=item *

absorbing_state

=item *

states_limit


=item *

info

=back


=head1 METHODS OF C<Bio::RNA::Treekin::PopulationDataRecord>

=head2 C<()>

TODO add


=head1 AUTHOR

Felix Kuehnl, C<< <felix@bioinf.uni-leipzig.de> >>


=head1 BUGS

Please report any bugs or feature requests to C<bug-bio-rna-treekin at
rt.cpan.org>, or through the web interface at
L<https://rt.cpan.org/NoAuth/ReportBug.html?Queue=Bio-RNA-Treekin>.  I will be
notified, and then you'll automatically be notified of progress on your bug as
I make changes.


=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc Bio::RNA::Treekin


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker (report bugs here)

L<https://rt.cpan.org/NoAuth/Bugs.html?Dist=Bio-RNA-Treekin>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/Bio-RNA-Treekin>

=item * CPAN Ratings

L<https://cpanratings.perl.org/d/Bio-RNA-Treekin>

=item * Search CPAN

L<https://metacpan.org/release/Bio-RNA-Treekin>

=back


# =head1 ACKNOWLEDGEMENTS


=head1 LICENSE AND COPYRIGHT

Copyright 2019 Felix Kuehnl.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see L<http://www.gnu.org/licenses/>.


=cut

1; # End of Bio::RNA::Treekin

