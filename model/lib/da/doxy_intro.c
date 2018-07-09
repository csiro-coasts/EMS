/**
 * \mainpage Data Assimilation library
 * \brief This DA library is a C++ implementaion of the Ensemble
 * Optimal Interpolation scheme of Jones et al [2012].
 *
 * \n \n
 * The library uses the GNU Scientific library (gsl) for all
 * matrix and vector manipulations, in particular, matrix inversion.
 * \n
 * 
 * The main output is:\n
 * \f[
 *    w^a=w^b + K\left(w^o - Hw^b\right)
 * \f]
 * The superscripts, \f$a\f$, \f$o\f$ and \f$b\f$ represent the
 * analysis, observation and background fields of the \f$w\f$ vector,
 * respectively, which is defined as:\n
 * \f[
 * w =
 * \left[\begin{array}{c}v1_1\\v1_2\\...\\v1_{n1}\\v2_1\\v2_2\\...\\v2_{n2}\\...\\vm_{nm}\end{array}\right]
 * \f]
 * where \f$v1...vm\f$ are the \f$m\f$ state variables that are being
 * assimilated and \f$n1...nm\f$ are the associated cell locations. Note,
 * typically these will only be either the set of all 3D or 2D wet cells
 *
 * \n
 * The Kalman gain is:
 * \n
 * \f[
 *   K = PH^t\left(HPH^t + R\right)^{-1}
 * \f]
 * where \f$H\f$ is an operator which maps locations in space from observation
 * to grid
 * \f[
 *   P = \frac{AA^t}{n-1}
 * \f]
 *
 * The C interface is availabe at da_interface.cpp.  Here is the basic
 * road map on how to use the library. The DA library has no concept
 * of cartesian space, it only deals with its own local sparse array
 * which is constructed by the caller as a stacked 1D array of
 * states. This stacked state array is what is used to construct the A
 * matrix as well as the background state vector \f$w^{b}\f$
 *
 * 1. Create the main object using da_create_object()<br>
 * 2. Parse the input prm-file to get sizes for each state and add using da_add_state()<br>
 * 3. Allocate the A matrix with da_allocA() <br>
 * 4. Construct local concept of valid sparse cells <br> 
 * 5. Fill in the A matrix using either da_fillA() or da_fillA_col()<br>
 * 6. Create maps using da_create_maps() and use da_fill_map() to fill
 * in 2D/3D cartesian index to local sparse maps<br>
 * 7. Initaialise observations from prm-file <br>
 * 8. Initialisation is now done <br>
 * 9. Begin sim-loop: <br>
 * 9a. Fill in background vector \f$w^b\f$ using
 * da_set_wb_col()/da_set_wb()
 * 9b. Call da_do_analysis() to perform the analysis
 * 9c. Retrieve \f$w^a\f$ using da_fill_wa()/da_get_wa()
 *
 * Step 9 can be repeated for different times.
 */
