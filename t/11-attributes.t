#!perl -T
use 5.012;
use warnings;

use Test::More;
use Test::NoWarnings;           # produces one additional test!
use File::Spec::Functions qw(catfile);
use Bio::RNA::Treekin;

plan tests => 5;


##############################################################################
##                                Input data                                ##
##############################################################################

my $treekin_single_small = catfile qw(t data treekin_single_small.kin);
my $treekin_multi_small  = catfile qw(t data treekin_multi_small.kin);

my $pops_single_small = [
    [qw( 1.00000000000000005551e-01 9.975614e-01 1.780033e-03 1.699285e-04 4.366637e-04 2.853783e-05 1.732678e-05 5.203420e-06 2.284726e-07 7.195039e-07 )],
    [qw( 3.29999999999999960032e-01 9.925545e-01 5.493689e-03 5.726556e-04 1.206952e-03 9.145588e-05 5.531990e-05 1.686918e-05 2.355482e-06 6.251624e-06 )],
    [qw( 8.58999999999999874767e-01 9.835173e-01 1.233819e-02 1.520550e-03 2.184748e-03 2.229225e-04 1.337191e-04 4.220111e-05 1.405266e-05 2.627692e-05 )],
    [qw( 2.07569999999999943441e+00 9.707329e-01 2.190646e-02 3.592723e-03 2.808839e-03 4.666975e-04 2.750026e-04 9.350034e-05 6.143765e-05 6.249221e-05 )],
    [qw( 4.87410999999999816623e+00 9.586998e-01 2.924071e-02 7.377603e-03 2.945105e-03 8.148683e-04 4.643097e-04 1.828670e-04 1.859488e-04 8.880556e-05 )],
    [qw( 1.13104529999999954271e+01 9.513184e-01 3.080657e-02 1.240492e-02 3.020511e-03 1.099316e-03 5.992051e-04 2.952576e-04 3.610896e-04 9.469268e-05 )],
    [qw( 2.61140418999999894822e+01 9.475712e-01 3.071786e-02 1.589999e-02 3.076343e-03 1.171927e-03 6.248600e-04 3.656116e-04 4.772013e-04 9.499057e-05 )],
    [qw( 6.01622963699999715459e+01 9.468848e-01 3.069401e-02 1.656599e-02 3.087073e-03 1.172795e-03 6.246081e-04 3.765550e-04 4.991709e-04 9.502509e-05 )],
    [qw( 1.38473281650999922476e+02 9.468742e-01 3.069365e-02 1.657630e-02 3.087240e-03 1.172769e-03 6.245952e-04 3.766629e-04 4.995109e-04 9.502564e-05 )],
    [qw( 1.00000000000000000000e+12 9.468742e-01 3.069365e-02 1.657630e-02 3.087240e-03 1.172769e-03 6.245952e-04 3.766629e-04 4.995109e-04 9.502564e-05 )],
];

my $pops_multi_small = [
    [qw( 1.00000000000000005551e-01 9.996821e-01 5.125677e-07 3.173154e-04 3.640100e-08 1.693872e-08 1.422256e-08 2.222878e-08 2.830810e-09 )],
    [qw( 2.36000000000000043077e-01 9.992505e-01 2.844770e-06 7.461457e-04 2.013262e-07 9.349781e-08 7.719580e-08 1.196145e-07 1.553593e-08 )],
    [qw( 4.20960000000000000853e-01 9.986650e-01 9.008105e-06 1.324396e-03 6.345097e-07 2.938785e-07 2.372402e-07 3.634494e-07 4.845561e-08 )],
    [qw( 6.72505600000000036687e-01 9.978714e-01 2.284224e-05 2.101812e-03 1.598693e-06 7.377558e-07 5.780181e-07 8.725282e-07 1.203874e-07 )],
    [qw( 1.01460761600000015648e+00 9.967971e-01 5.154067e-05 3.142840e-03 3.576117e-06 1.642208e-06 1.236886e-06 1.832157e-06 2.642883e-07 )],
    [qw( 1.47986635776000019504e+00 9.953448e-01 1.083623e-04 4.529507e-03 7.430967e-06 3.390019e-06 2.425411e-06 3.508946e-06 5.356373e-07 )],
    [qw( 2.11261824655360053171e+00 9.933860e-01 2.173581e-04 6.363655e-03 1.467111e-05 6.634396e-06 4.443921e-06 6.248833e-06 1.023265e-06 )],
    [qw( 2.97316081531289677642e+00 9.907510e-01 4.214040e-04 8.767371e-03 2.784320e-05 1.244500e-05 7.674853e-06 1.044426e-05 1.860360e-06 )],
    [qw( 4.14349870882553972251e+00 9.872192e-01 7.953769e-04 1.187968e-02 5.106941e-05 2.247964e-05 1.253746e-05 1.646944e-05 3.229337e-06 )],
    [qw( 5.73515824400273466210e+00 9.825072e-01 1.466712e-03 1.584693e-02 9.064509e-05 3.911707e-05 1.940148e-05 2.461175e-05 5.349154e-06 )],
    [qw( 7.89981521184371970890e+00 9.762576e-01 2.645159e-03 2.080442e-02 1.554147e-04 6.540187e-05 2.849307e-05 3.504259e-05 8.432535e-06 )],
    [qw( 1.08437486881074605094e+01 9.680288e-01 4.660932e-03 2.684912e-02 2.563175e-04 1.045690e-04 3.985043e-05 4.782904e-05 1.261317e-05 )],
    [qw( 1.48474982158261461507e+01 9.572896e-01 8.003640e-03 3.400946e-02 4.042276e-04 1.589689e-04 5.334690e-05 6.292888e-05 1.787440e-05 )],
    [qw( 2.02925975735235581965e+01 9.434222e-01 1.334185e-02 4.222883e-02 6.055982e-04 2.286315e-04 6.873406e-05 8.015466e-05 2.403838e-05 )],
    [qw( 2.76979326999920409946e+01 9.257404e-01 2.148825e-02 5.138774e-02 8.573436e-04 3.104967e-04 8.569174e-05 9.919701e-05 3.084884e-05 )],
    [qw( 3.77691884719891817213e+01 9.035397e-01 3.327168e-02 6.138155e-02 1.145519e-03 3.996818e-04 1.039555e-04 1.197849e-04 3.810462e-05 )],
    [qw( 5.14660963219052902673e+01 8.762111e-01 4.930797e-02 7.222452e-02 1.452582e-03 4.926591e-04 1.235135e-04 1.419174e-04 4.576289e-05 )],
    [qw( 7.00938909977912061322e+01 8.434695e-01 6.972998e-02 8.407673e-02 1.769725e-03 5.895118e-04 1.446696e-04 1.659330e-04 5.394728e-05 )],
    [qw( 9.54276917569960403398e+01 8.057303e-01 9.398183e-02 9.707278e-02 2.099986e-03 6.923689e-04 1.677275e-04 1.921547e-04 6.280424e-05 )],
    [qw( 1.29881660789514626231e+02 7.645669e-01 1.207194e-01 1.109830e-01 2.444852e-03 8.010456e-04 1.923491e-04 2.201745e-04 7.223563e-05 )],
    [qw( 1.76739058673739918959e+02 7.230028e-01 1.477860e-01 1.249649e-01 2.789239e-03 9.099364e-04 2.170834e-04 2.483273e-04 8.170375e-05 )],
    [qw( 2.40465119796286302289e+02 6.852424e-01 1.723842e-01 1.376593e-01 3.101627e-03 1.008757e-03 2.395383e-04 2.738863e-04 9.029855e-05 )],
    [qw( 3.27132562922949375661e+02 6.555128e-01 1.917514e-01 1.476535e-01 3.347552e-03 1.086556e-03 2.572167e-04 2.940086e-04 9.706507e-05 )],
    [qw( 4.45000285575211194100e+02 6.362148e-01 2.043229e-01 1.541408e-01 3.507184e-03 1.137056e-03 2.686920e-04 3.070703e-04 1.014573e-04 )],
    [qw( 6.05300388382287223976e+02 6.265511e-01 2.106183e-01 1.573895e-01 3.587123e-03 1.162344e-03 2.744385e-04 3.136111e-04 1.036568e-04 )],
    [qw( 8.23308528199910711010e+02 6.231324e-01 2.128453e-01 1.585387e-01 3.615402e-03 1.171291e-03 2.764714e-04 3.159250e-04 1.044349e-04 )],
    [qw( 1.00000000000000000000e+03 6.225084e-01 2.132519e-01 1.587485e-01 3.620565e-03 1.172924e-03 2.768425e-04 3.163474e-04 1.045769e-04 )],
];

##############################################################################
##                              Test functions                              ##
##############################################################################

# Given an object and a key--value list of attributes and values, check that
# the object has a method called <attribute> and that it returns <value>.
sub test_attrib_vals {
    my $obj           = shift;
    my $test_name     = shift;
    my @val_of_attrib = @_;

    subtest $test_name => sub {
        plan tests => 2 * @val_of_attrib/2;
        while (my ($attrib, $val) = splice @val_of_attrib, 0, 2) {
            can_ok $obj, $attrib;
            is $obj->$attrib, $val, "val of $attrib";
        }
    };

}

sub test_pop_data {
    my ($record, $record_name, $pops_data) = @_;

    my $population_data_count = @$pops_data;
    my $min_count             = @{ $pops_data->[0] } - 1;    # minus time field
    subtest "$record_name population data" => sub {
        plan tests => 2  +  $population_data_count * ($min_count + 1);

        is $record->population_data_count, $population_data_count, 'population_data_count';
        is $record->min_count, $min_count, 'min_count';

        # Iterate over population data rows.
        foreach my $i (0..($population_data_count-1)) {
            my $pop          = $record->population($i);
            my $pop_expected = $pops_data->[$i];
            is $pop->time, $pop_expected->[0], "data row $i: times";

            # Iterate over population values.
            foreach my $min (1..$min_count) {
                is $pop->of_min($min), $pop_expected->[$min], "data row $i: pop of min $min";
            }
        }
    };
}

##############################################################################
##                                Call tests                                ##
##############################################################################


my $record_single_small
    = Bio::RNA::Treekin::Record->new($treekin_single_small);
test_attrib_vals $record_single_small,
                 'record attribs single small',
                 'date'              => 'Wed Feb 12 11:23:22 2020',
                 'sequence'          => '(null)',
                 'method'            => 'I',
                 'start_time'        => '0.10',
                 'stop_time'         => '1000000000000.00',
                 'temperature'       => '37.00',
                 'basename'          => '<stdin>',
                 'time_increment'    => '2.30',
                 'degeneracy'        => ' off',
                 'absorbing_state'   => 'none',
                 'states_limit'      => '2147483647',
                 ;
test_pop_data $record_single_small, 'single small', $pops_single_small;

# Get only first record of multi-file
my $record_multi_small
    = Bio::RNA::Treekin::MultiRecord->new($treekin_multi_small)->next();
test_attrib_vals $record_multi_small,
                 'record attribs multi small',
                 'date'             => 'Wed Feb 12 11:52:49 2020',
                 'sequence'         => '(null)',
                 'method'           => 'I',
                 'start_time'       => '0.10',
                 'stop_time'        => '1000.00',
                 'temperature'      => '37.00',
                 'basename'         => '<stdin>',
                 'time_increment'   => '1.36',
                 'degeneracy'       => ' off',
                 'absorbing_state'  => 'none',
                 'states_limit'     => '2147483647',
                 # Multi-exclusive attributes:
                 'rates_file'       => '28.rates.bin',
                 'file_index'       => '0',
                 'cmd'              => 'treekin -m I --feps=-1.0 --tinc=1.36 --t0=0.1 --t8=1e3 --Temp=37 --bin --recoverE < 28.rates.bin',
                 ;
test_pop_data $record_multi_small, 'multi small', $pops_multi_small;

exit 0;                             # EOF