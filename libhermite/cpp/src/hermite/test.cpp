// #include <boost/python.hpp>
// #include <boost/python/numpy.hpp>
// #include <boost/python/suite/indexing/vector_indexing_suite.hpp>

// using namespace std;

// typedef vector<double> vec;
// typedef vector<vec> mat;

// mat test()
// {
//     int n = 5e3;
//     mat result(n, vec(n, 0.));
//     for (u_int i = 0; i < result.size(); i++)
//     {
//         for (u_int j = 0; j < result[0].size(); j++)
//         {
//             result[i][j] = (double) i - (double) j;
//         }
//     }
//     return result;
// }

// namespace p = boost::python;
// namespace np = boost::python::numpy;

// np::ndarray convert_to_numpy()
// {
//     // u_int n_rows = input.size();
//     // u_int n_cols = input[0].size();
//     // p::tuple shape = p::make_tuple(n_rows, n_cols);
//     // np::dtype dtype = np::dtype::get_builtin<float>();
//     // p::tuple stride = p::make_tuple(n_cols,1);
//     // p::object own;
//     // np::ndarray converted = np::from_data(input.data(), dtype, shape, stride, own);
//     int data[] = {1,2,3,4,5};
//     p::tuple shape = p::make_tuple(5);
//     p::tuple stride = p::make_tuple(sizeof(int));
//     p::object own;
//     np::dtype dt = np::dtype::get_builtin<int>();
//     np::ndarray converted = np::from_data(data, dt, shape,stride,own);
//     return converted;
// }


// BOOST_PYTHON_MODULE(hermite_cpp)
// {
//     using namespace boost::python;

//     // Initialize numpy
//     Py_Initialize();
//     boost::python::numpy::initialize();

//     class_<vec>("double_vec")
//         .def(vector_indexing_suite<vec>())
//         ;

//     class_<mat>("double_mat")
//         .def(vector_indexing_suite<mat>())
//         ;

//     def("convert_to_numpy", convert_to_numpy);
//     def("test", test);
// }
