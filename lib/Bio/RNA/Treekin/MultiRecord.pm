# Bio/RNA/Treekin/MultiRecord.pm
package Bio::RNA::Treekin::MultiRecord;
our $VERSION = '0.01';

use 5.006;
use strict;
use warnings;

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

1; # End of Bio::RNA::Treekin::MultiRecord


__END__


=pod

=encoding UTF-8

=head1 NAME

Bio::RNA::Treekin::MultiRecord - TODO I<Treekin> output.

=head1 SYNOPSIS

    use Bio::RNA::Treekin;

=head1 DESCRIPTION

=head1 METHODS

=head2 new()

=head2 ...()


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

# End of Bio/RNA/Treekin/MultiRecord.pm
