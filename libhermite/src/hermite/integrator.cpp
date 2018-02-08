#include <cmath>
#include <functional>
#include <numeric>
#include <iostream>
#include <vector>

#include <boost/function.hpp>
#include <boost/python.hpp>
#include <boost/math/special_functions/binomial.hpp>

#include "hermite/helpers/templates.hpp"
#include "hermite/helpers/combinatorics.hpp"
#include "hermite/integrator.hpp"

#define PI 3.141592653589793238462643383279502884

typedef std::vector< std::vector< std::vector<double> > > cube;
typedef std::vector< std::vector<double> > mat;
typedef std::vector<double> vec;

using namespace std;

namespace HermiteQuadrature {

    static vector<double> nodes_6 = { 2.35060497367449222281e+00, 4.36077411927616508688e-01, 1.33584907401369694976e+00 };
    static vector<double> weights_6 = { 4.53000990550884564102e-03, 7.24629595224392524086e-01, 1.57067320322856643914e-01 };

    static vector<double> nodes_10 = { 3.42901327223704608781e-01, 1.03661082978951365422e+00, 1.75668364929988177344e+00, 2.53273167423278979644e+00, 3.43615911883773760341e+00 };
    static vector<double> weights_10 = { 6.10862633735325798764e-01, 2.40138611082314686412e-01, 3.38743944554810631376e-02, 1.34364574678123269223e-03, 7.64043285523262062930e-06 };

    static vector<double> nodes_20 = { 2.45340708300901249903e-01, 7.37473728545394358719e-01, 1.23407621539532300786e+00, 1.73853771211658620678e+00, 2.25497400208927552311e+00, 2.78880605842813048055e+00, 3.34785456738321632688e+00, 3.94476404011562521040e+00, 4.60368244955074427298e+00, 5.38748089001123286199e+00 };
    static vector<double> weights_20 = { 4.62243669600610089640e-01, 2.86675505362834129720e-01, 1.09017206020023320014e-01, 2.48105208874636108814e-02, 3.24377334223786183217e-03, 2.28338636016353967260e-04, 7.80255647853206369398e-06, 1.08606937076928169398e-07, 4.39934099227318055366e-10, 2.22939364553415129254e-13 };

    static vector<double> nodes_30 = { 0.2011285765488714855458, 0.6039210586255523077782, 1.008338271046723461805, 1.415527800198188511941, 1.826741143603688038836, 2.243391467761504072473, 2.667132124535617200571, 3.099970529586441748689, 3.544443873155349886925, 4.003908603861228815228, 4.483055357092518341887, 4.988918968589943944486, 5.533147151567495725118, 6.138279220123934620395, 6.863345293529891581061};
    static vector<double> weights_30 = { 0.3863948895418138625556, 0.2801309308392126674135, 0.1467358475408900997517, 0.05514417687023425116808, 0.01470382970482668351528, 0.00273792247306765846299, 3.48310124318685523421E-4, 2.9387252289229876415E-5, 1.57909488732471028835E-6, 5.10852245077594627739E-8, 9.178580424378528209E-10, 8.10618629746304420399E-12, 2.87860708054870606219E-14, 2.8103336027509037088E-17, 2.90825470013122622941E-21};

    static vector<double> nodes_100 = { 1.10795872422439482889e-01, 3.32414692342231807054e-01, 5.54114823591616988249e-01, 7.75950761540145781976e-01, 9.97977436098105243902e-01, 1.22025039121895305882e+00, 1.44282597021593278768e+00, 1.66576150874150946983e+00, 1.88911553742700837153e+00, 2.11294799637118795206e+00, 2.33732046390687850509e+00, 2.56229640237260802502e+00, 2.78794142398198931316e+00, 3.01432358033115551667e+00, 3.24151367963101295043e+00, 3.46958563641858916968e+00, 3.69861685931849193984e+00, 3.92868868342767097205e+00, 4.15988685513103054019e+00, 4.39230207868268401677e+00, 4.62603063578715577309e+00, 4.86117509179121020995e+00, 5.09784510508913624692e+00, 5.33615836013836049734e+00, 5.57624164932992410311e+00, 5.81823213520351704715e+00, 6.06227883261430263882e+00, 6.30854436111213512156e+00, 6.55720703192153931598e+00, 6.80846335285879641431e+00, 7.06253106024886543766e+00, 7.31965282230453531632e+00, 7.58010080785748888415e+00, 7.84418238446082116862e+00, 8.11224731116279191689e+00, 8.38469694041626507474e+00, 8.66199616813451771409e+00, 8.94468921732547447845e+00, 9.23342089021916155069e+00, 9.52896582339011480496e+00, 9.83226980777796909401e+00, 1.01445099412928454695e+01, 1.04671854213428121416e+01, 1.08022607536847145950e+01, 1.11524043855851252649e+01, 1.15214154007870302416e+01, 1.19150619431141658018e+01, 1.23429642228596742953e+01, 1.28237997494878089065e+01, 1.34064873381449101387e+01 };
    static vector<double> weights_100 = { 2.18892629587439125060e-01, 1.98462850254186477710e-01, 1.63130030502782941425e-01, 1.21537986844104181985e-01, 8.20518273912244646789e-02, 5.01758126774286956964e-02, 2.77791273859335142698e-02, 1.39156652202318064178e-02, 6.30300028560805254921e-03, 2.57927326005909017346e-03, 9.52692188548619117497e-04, 3.17291971043300305539e-04, 9.51716277855096647040e-05, 2.56761593845490630553e-05, 6.22152481777786331722e-06, 1.35179715911036728661e-06, 2.62909748375372507934e-07, 4.56812750848493951350e-08, 7.07585728388957290740e-09, 9.74792125387162124528e-10, 1.19130063492907294976e-10, 1.28790382573155823282e-11, 1.22787851441012497000e-12, 1.02887493735099254677e-13, 7.54889687791524329227e-15, 4.82983532170303334787e-16, 2.68249216476037608006e-17, 1.28683292112115327575e-18, 5.30231618313184868536e-20, 1.86499767513025225814e-21, 5.56102696165916731717e-23, 1.39484152606876708047e-24, 2.91735007262933241788e-26, 5.03779116621318778423e-28, 7.10181222638493422964e-30, 8.06743427870937717382e-32, 7.27457259688776757460e-34, 5.11623260438522218054e-36, 2.74878488435711249209e-38, 1.10047068271422366943e-40, 3.18521787783591793076e-43, 6.42072520534847248278e-46, 8.59756395482527161007e-49, 7.19152946346337102982e-52, 3.45947793647555044453e-55, 8.51888308176163378638e-59, 9.01922230369355617950e-63, 3.08302899000327481204e-67, 1.97286057487945255443e-72, 5.90806786503120681541e-79 };

    static vector<double> nodes_2 = {7.07106781186547524382e-01};
    static vector<double> weights_2 = {8.86226925452758013655e-01};

    static vector<double> nodes_4 = { 5.24647623275290317900e-01, 1.65068012388578455585e+00 };
    static vector<double> weights_4 = { 8.04914090005512836482e-01, 8.13128354472451771398e-02 };

    static vector<double> nodes_8 = { 3.81186990207322116844e-01, 1.15719371244678019474e+00, 1.98165675669584292584e+00, 2.93063742025724401920e+00 };
    static vector<double> weights_8 = { 6.61147012558241291042e-01, 2.07802325814891879546e-01, 1.70779830074134754563e-02, 1.99604072211367619211e-04 };

    static vector<double> nodes_16 = { 2.73481046138152452172e-01, 8.22951449144655892596e-01, 1.38025853919888079639e+00, 1.95178799091625397740e+00, 2.54620215784748136221e+00, 3.17699916197995602682e+00, 3.86944790486012269869e+00, 4.68873893930581836465e+00 };
    static vector<double> weights_16 = { 5.07929479016613741923e-01, 2.80647458528533675357e-01, 8.38100413989858294132e-02, 1.28803115355099736832e-02, 9.32284008624180529895e-04, 2.71186009253788151199e-05, 2.32098084486521065344e-07, 2.65480747401118224476e-10 };

    static vector<double> nodes_32 = { 1.94840741569399326713e-01, 5.84978765435932448449e-01, 9.76500463589682838499e-01, 1.37037641095287183817e+00, 1.76765410946320160465e+00, 2.16949918360611217335e+00, 2.57724953773231745414e+00, 2.99249082500237420621e+00, 3.41716749281857073593e+00, 3.85375548547144464390e+00, 4.30554795335119844506e+00, 4.77716450350259639289e+00, 5.27555098651588012760e+00, 5.81222594951591383294e+00, 6.40949814926966041214e+00, 7.12581390983072757292e+00 };
    static vector<double> weights_32 = { 3.75238352592802392864e-01, 2.77458142302529898131e-01, 1.51269734076642482578e-01, 6.04581309559126141860e-02, 1.75534288315734303030e-02, 3.65489032665442807915e-03, 5.36268365527972045989e-04, 5.41658406181998255789e-05, 3.65058512956237605727e-06, 1.57416779254559402923e-07, 4.09883216477089661816e-09, 5.93329146339663861478e-11, 4.21501021132644757306e-13, 1.19734401709284866582e-15, 9.23173653651829223381e-19, 7.31067642738416239302e-23 };

    static vector<double> nodes_64 = { 1.38302244987009724116e-01, 4.14988824121078684584e-01, 6.91922305810044577278e-01, 9.69269423071178016745e-01, 1.24720015694311794072e+00, 1.52588914020986366295e+00, 1.80551717146554491895e+00, 2.08627287988176202084e+00, 2.36835458863240140418e+00, 2.65197243543063501106e+00, 2.93735082300462180976e+00, 3.22473129199203572583e+00, 3.51437593574090621144e+00, 3.80657151394536046127e+00, 4.10163447456665671486e+00, 4.39991716822813764762e+00, 4.70181564740749981593e+00, 5.00777960219876819639e+00, 5.31832522463327085713e+00, 5.63405216434997214731e+00, 5.95566632679948604551e+00, 6.28401122877482823555e+00, 6.62011226263602737922e+00, 6.96524112055110752916e+00, 7.32101303278094920138e+00, 7.68954016404049682824e+00, 8.07368728501022522599e+00, 8.47752908337986309066e+00, 8.90724909996476975696e+00, 9.37315954964672116284e+00, 9.89528758682953902149e+00, 1.05261231679605458834e+01 };
    static vector<double> weights_64 = { 2.71377424941303977947e-01, 2.32994786062678046649e-01, 1.71685842349083702001e-01, 1.08498349306186840632e-01, 5.87399819640994345496e-02, 2.72031289536889184544e-02, 1.07560405098791370492e-02, 3.62258697853445876057e-03, 1.03632909950757766343e-03, 2.50983698513062486093e-04, 5.12592913578627466080e-05, 8.78849923085035918142e-06, 1.25834025103118457615e-06, 1.49553293672724706114e-07, 1.46512531647610935494e-08, 1.17361674232154934349e-09, 7.61521725014545135318e-11, 3.95917776694772392711e-12, 1.62834073070972036212e-13, 5.21862372659084752279e-15, 1.28009339132243804163e-16, 2.35188471067581911698e-18, 3.15225456650378141599e-20, 2.98286278427985115448e-22, 1.91170688330064282993e-24, 7.86179778892591036891e-27, 1.92910359546496685035e-29, 2.54966089911299925654e-32, 1.55739062462976380226e-35, 3.42113801125574050436e-39, 1.67974799010815921869e-43, 5.53570653585694282057e-49 };

    // nodes and weight shuold have the right length
    void get_gh_quadrature(int nNodes, vector<double>& nodes, vector<double>& weights) {

        nodes = vec(nNodes);
        weights = vec(nNodes);

        // Unscaled nodes and weights
        vector<double> x;
        vector<double> w;

        switch (nNodes) {
            case 1: nodes = {0.}; weights = {1.}; return; break;
            case 2: x = nodes_2; w = weights_2; break;
            case 4: x = nodes_4; w = weights_4; break;
            case 6: x = nodes_6; w = weights_6; break;
            case 8: x = nodes_8; w = weights_8; break;
            case 10: x = nodes_10; w = weights_10; break;
            case 16: x = nodes_16; w = weights_16; break;
            case 20: x = nodes_20; w = weights_20; break;
            case 30: x = nodes_30; w = weights_30; break;
            case 32: x = nodes_32; w = weights_32; break;
            case 64: x = nodes_64; w = weights_64; break;
            case 100: x = nodes_100; w = weights_100; break;
            default: cout << "Invalid number of nodes (" << nNodes << ") for Gauss-hermite integration" << endl; exit(0);
        }

        // Add symmetric
        for (int i = 0; i < nNodes/2; ++i) {
            nodes[i] = x[i];
            weights[i] = w[i];
            nodes[i + nNodes/2] = -x[i];
            weights[i + nNodes/2] = w[i];
        }

        // Scaling to get rid of factors
        for (int i = 0; i < nNodes; ++i) {
            weights[i] /= sqrt(PI);
            nodes[i] *= sqrt(2);
        }
    }

    // Compute the nodes and weights resulting from the product of quadratures
    void quad_prod(vector<int> sizes, mat& nodes, vec& weights) {

        // Dimension
        int dim = sizes.size();

        // Number of nodes of product
        int n = 1;
        for (int i = 0; i < dim; ++i)
            n *= sizes[i];

        // Store nodes and weights of all the quadrature rules used
        mat all_nodes(dim);
        mat all_weights(dim);

        for (int i = 0; i < dim; ++i)
            get_gh_quadrature(sizes[i], all_nodes[i], all_weights[i]);

        // Initialize
        nodes = mat (n, vec(dim, 0.));
        weights = vec (n, 1.);

        // Compute nodes and weights of product
        for (int i = 0; i < n; ++i) {
            for (int j = 0, aux = i; j < dim; ++j) {
                nodes[i][j] = all_nodes[j][aux % sizes[j]];
                weights[i] *= all_weights[j][aux % sizes[j]];
                aux = aux / sizes[j];
            }
        }
    }

    double integrate_nd(int nNodes, int dim, boost::function<double(vec)> func) {
        mat nodes;
        vec weights;

        vector<int> sizes(dim, nNodes);
        quad_prod(sizes, nodes, weights);

        double result = 0.;
        for (unsigned int i = 0; i < nodes.size(); ++i) {
            result += func(nodes[i]) * weights[i];
        }
        return result;
    }

    double integrate(int nNodes, boost::function<double(double)> func) {

        vec nodes;
        vec weights;
        get_gh_quadrature(nNodes, nodes, weights);

        double result = 0.;
        for (unsigned int i = 0; i < nodes.size(); ++i) {
            result += func(nodes[i]) * weights[i];
        }
        return result;
    }

    void get_smolyak_nodes_and_weights(int nVars, vector< vector<double> >& snodes, vector<double>& sweights) {

        snodes = vector< vector<double> > (0);
        sweights = vector<double>(0);

        vector<int> indices = {0, 1, 2, 3, 4, 5, 6};
        vector<int> sizes = {1, 2, 4, 8, 16, 32, 64};

        // Using notations of
        // https://dbarajassolano.wordpress.com/2012/01/26/on-sparse-grid-quadratures/
        int q = indices.size() - 1;

        // Enumeration of multi_indinces
        vector< vector<int> > enum_indices = interval_multi_indices(nVars, q - nVars + 1, q);

        // Auxiliary vectors
        vector< vector<double> > aux_nodes;
        vector<double> aux_weights;

        for (unsigned int i = 0; i < enum_indices.size(); ++i) {

            vector<int> this_index = enum_indices[i];
            vector<int> this_sizes(this_index.size());

            for (unsigned int j = 0; j < this_index.size(); ++j)
                this_sizes[j] = sizes[this_index[j]];

            int sum_index = accumulate(this_index.begin(), this_index.end(), 0.);
            double weight_mod = pow(-1, q - sum_index) * boost::math::binomial_coefficient<double> (nVars - 1, q - sum_index);

            quad_prod(this_sizes, aux_nodes, aux_weights);
            aux_weights = aux_weights * weight_mod;

            snodes.insert(snodes.end(), aux_nodes.begin(), aux_nodes.end());
            sweights.insert(sweights.end(), aux_weights.begin(), aux_weights.end());
        }
    }

    class Quad {
        public: 
            Quad(int, int);
            mat nodes;
            vec weights;
    };

    HermiteQuadrature::Quad::Quad(int nNodes, int nVars) {

        vector<double> nodes_1d(nNodes);
        vector<double> weights_1d(nNodes);

        if (nNodes == 0)
            get_smolyak_nodes_and_weights(nVars, nodes, weights);
        else {
            vector<int> seq(nVars, nNodes);
            quad_prod(seq, nodes, weights);
        }
    }

    // ---- PYTHON WRAPPERS ----
    double integrate_wrap(int nNodes, boost::python::object func)
    {
        auto lambda = [func](double a) -> double {
            return boost::python::extract<double>(func(a));
        };
        return integrate(nNodes, boost::function<double(double)>(lambda));
    }

    /*!
     * Integrate 2D function
     */
    double integrate_2d_wrap(int nNodes, boost::python::object func)
    {
        auto lambda = [func](vec v) -> double {
            return boost::python::extract<double>(func(v[0],v[1]));
        };
        return integrate_nd(nNodes, 2, boost::function<double(vec)>(lambda));
    }

    double integrate_3d_wrap(int nNodes, boost::python::object func)
    {
        auto lambda = [func](vec v) -> double {
            return boost::python::extract<double>(func(v[0],v[1],v[2]));
        };
        return integrate_nd(nNodes, 3, boost::function<double(vec)>(lambda));
    }

    // ---- EXPOSE TO PYTHON ----
    BOOST_PYTHON_MODULE(hermite)
    {
        using namespace boost::python;
        def("integrate", integrate_wrap);
        def("integrate_2d", integrate_2d_wrap);
        def("integrate_3d", integrate_3d_wrap);
    }
}
