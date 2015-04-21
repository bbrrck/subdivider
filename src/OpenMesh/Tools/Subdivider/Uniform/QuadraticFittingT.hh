/*                                                                           *
 *                               OpenMesh                                    *
 *      Copyright (C) 2001-2013 by Computer Graphics Group, RWTH Aachen      *
 *                           www.openmesh.org                                *
 *                                                                           */
//=============================================================================
//
//  CLASS QuadraticFittingT
//
//		OpenMesh implementation of Quadric Fitting Refinement :
//		nonlinear interpolating subdivision scheme for triangle mesh
//
//                 author:   Tibor Stanko
//		last modification:   April 01, 2014, 17:35
//
//=============================================================================

#ifndef OPENMESH_SUBDIVIDER_UNIFORM_QUADRATICFITTINGT_HH
#define OPENMESH_SUBDIVIDER_UNIFORM_QUADRATICFITTINGT_HH

//== INCLUDES =================================================================

#include <OpenMesh/Core/System/config.hh>
#include <OpenMesh/Tools/Subdivider/Uniform/SubdividerT.hh>
#include <OpenMesh/Core/Utils/vector_cast.hh>
// -------------------- STL
#include <vector>
#include <algorithm>
#if defined(OM_CC_MIPS)
#  include <math.h>
#else
#  include <cmath>
#endif
// #include <unordered_set>
#include <armadillo>

//== NAMESPACE ================================================================

namespace OpenMesh   { // BEGIN_NS_OPENMESH
namespace Subdivider { // BEGIN_NS_SUBDIVIDER
namespace Uniform    { // BEGIN_NS_UNIFORM


//== CLASS DEFINITION =========================================================

template <typename MeshType, typename RealType = float>
class QuadraticFittingT : public SubdividerT<MeshType, RealType>
{
public:

	typedef RealType										real_t;
	typedef MeshType										mesh_t;
	typedef SubdividerT< mesh_t, real_t >					parent_t;

public:

	QuadraticFittingT(void) : parent_t() {}
	QuadraticFittingT( mesh_t& _m ) : parent_t(_m) {}
	~QuadraticFittingT() {}

	const char *name() const { return "Quadratic Fitting"; }
	void update_normals( mesh_t& _m ) {_m.update_face_normals();}

protected:

/*
	-----------
	PREPARATION
	-----------
*/
	bool prepare( mesh_t& _m )
	{
		_m.request_edge_status();
		_m.request_vertex_status();
		_m.add_property( vertex_level_ );
		_m.add_property( new_vertex_pos_ );
		_m.add_property( new_vertex_normal_ );
		_m.add_property( quadric_coeff_ );

		_m.request_vertex_normals();
		if( !_m.has_vertex_normals() ) 
		{
			_m.request_face_normals();
			_m.update_normals();
		}
		return _m.has_edge_status() && _m.has_vertex_normals() && _m.has_vertex_status() && vertex_level_.is_valid() 
				&& new_vertex_pos_.is_valid() && new_vertex_normal_.is_valid() && quadric_coeff_.is_valid();
	}


/*
	-------
	CLEANUP
	-------
*/
	bool cleanup( mesh_t& _m )
	{
		_m.release_edge_status();
		_m.release_vertex_normals();
		_m.release_face_normals();
		_m.release_vertex_status();
		_m.remove_property( vertex_level_ );
		_m.remove_property( new_vertex_pos_ );
		_m.remove_property( new_vertex_normal_ );
		_m.remove_property( quadric_coeff_ );
		return true;
	}


/*
	-----------
	SUBDIVISION
	-----------
*/
	bool subdivide( mesh_t& _m, size_t _n, const bool _update_points = true)
	{
		typedef std::set<typename mesh_t::VertexHandle> vhandles_t;

		typename mesh_t::FaceIter 		f_it, f_end = _m.faces_end();
		typename mesh_t::EdgeIter 		e_it;
		typename mesh_t::VertexIter 	v_it;
		typename MeshType::VertexHandle	vh;
	
		for (size_t i=0; i < _n; ++i) 
		{
			// output the used weights
			std::cout << "QFR!" << std::endl;
			std::cout << "weights:" << std::endl;
			std::cout << "\t vertex init = " << VERTEX_WEIGHT_INIT << std::endl;
			std::cout << "\t vertex fact = " << VERTEX_WEIGHT_FACTOR << std::endl;
			std::cout << "\t normal init = " << NORMAL_WEIGHT_INIT << std::endl;
			std::cout << "\t normal fact = " << NORMAL_WEIGHT_FACTOR << std::endl;
			std::cout << "\t normalize normals? " << NORMALIZE_NORMALS << std::endl;
			
			// mark old edges
			for (e_it=_m.edges_begin(); e_it != _m.edges_end();++e_it)
				_m.status( *e_it ).set_tagged( true );

			// compute positions of new vertices
			for( f_it = _m.faces_begin(); f_it != f_end; ++f_it )
				compute_new_vertex( _m, *f_it );

			// add new vertices and perform face split
			for( f_it = _m.faces_begin(); f_it != f_end; ++f_it )
			{
				vh = _m.add_vertex( _m.property( new_vertex_pos_, *f_it ) );
				_m.split( *f_it, vh );
				_m.set_normal( vh, _m.property( new_vertex_normal_, *f_it ) );
			}

			// flip old edges
			for ( e_it = _m.edges_begin(); e_it != _m.edges_end(); ++e_it )
				if ( _m.status( *e_it ).tagged() && !_m.is_boundary( *e_it ) && _m.is_flip_ok( *e_it ) )
					_m.flip( *e_it );
		}

		return true;
	}

private:

/*
	--------
	SETTINGS
	--------
*/
	// minimal cardinality for triangle neighborhood
	static const int MIN_CARDINALITY = 9;
	static const double CONTINUE_COUNT_MAX = 100;

	// Foot point algorithm (Hartmann)
	static const double EPSILON = 0.001;
	static const double ALPHA_MAX = 0.5;

	bool NORMALIZE_NORMALS;

	int VERTEX_PICK_METHOD;

	// weights
	float VERTEX_WEIGHT_INIT,
		  VERTEX_WEIGHT_FACTOR,
		  NORMAL_WEIGHT_INIT,
		  NORMAL_WEIGHT_FACTOR;
	/*
	classic:
		VERTEX_WEIGHT_INIT = 1		VERTEX_WEIGHT_FACTOR = 0.1		NORMAL_WEIGHT_INIT =  0.001		NORMAL_WEIGHT_FACTOR = 0.01

	spike sphere from cube:
		VERTEX_WEIGHT_INIT = 1		VERTEX_WEIGHT_FACTOR = 1		NORMAL_WEIGHT_INIT =  1			NORMAL_WEIGHT_FACTOR = 1

	ellipsoid from cube, quadric reconstruction:
		VERTEX_WEIGHT_INIT = 1000	VERTEX_WEIGHT_FACTOR = 1		NORMAL_WEIGHT_INIT =  0.001		NORMAL_WEIGHT_FACTOR = 1
	*/


public:
	
	void ChangeWeights(float new_vWi, float new_vWf, float new_nWi, float new_nWf, bool new_NN = true, int vertexPick = 2)
	{
		VERTEX_PICK_METHOD   = vertexPick;
		VERTEX_WEIGHT_INIT   = new_vWi;
		VERTEX_WEIGHT_FACTOR = new_vWf;
		NORMAL_WEIGHT_INIT   = new_nWi;// * 1000000;
		NORMAL_WEIGHT_FACTOR = new_nWf;
		NORMALIZE_NORMALS    = new_NN;
	}


private:

/*
	-----------
	HELPER DATA
	-----------
*/
	// distance of vertex from triangle
	OpenMesh::VPropHandleT< int > vertex_level_;

	// coordinates of the new vertex
	OpenMesh::FPropHandleT< typename mesh_t::Point > new_vertex_pos_;

	// normal at the new vertex
	OpenMesh::FPropHandleT< typename mesh_t::Point > new_vertex_normal_;

	// quadric - computed coefficients
	OpenMesh::FPropHandleT< arma::vec > quadric_coeff_;


/*
	---------------------
	TRIANGLE NEIGHBORHOOD
	---------------------
*/
	std::set<typename mesh_t::VertexHandle>* face_neighbors(mesh_t& _m, typename mesh_t::FaceHandle face_handle)
	{
		typedef std::set<typename mesh_t::VertexHandle> vhandles_t;
		typename vhandles_t::iterator v_it, v_end;
		typename mesh_t::VertexVertexIter vv_it;
		typename mesh_t::FaceVertexIter fv_it;

		vhandles_t* face_neighbors = new vhandles_t();
		int level = 1;

		// first, insert all face vertices to the list
		for(fv_it = _m.fv_iter(face_handle); fv_it.is_valid(); ++fv_it)
		{
			_m.property( vertex_level_, *fv_it ) = level;
			face_neighbors->insert( *fv_it );
		}

		int vertex_count = face_neighbors->size();
		int last_step = vertex_count;
		int continue_count = 0;
		v_end = face_neighbors->end();

		while( vertex_count < MIN_CARDINALITY ) 
		{
			level++;
			last_step = 0;

			for(v_it = face_neighbors->begin(); v_it != face_neighbors->end(); ++v_it)
			{
				if( _m.property( vertex_level_, *v_it) != level - 1 ) continue;
				for (vv_it = _m.vv_iter(*v_it); vv_it.is_valid(); ++vv_it)
				{
					if(!_m.status( *vv_it ).tagged())
						face_neighbors->insert( *vv_it );

					if( face_neighbors->size() > vertex_count ) // additional vertex was inserted, change its level
					{
						_m.property( vertex_level_, *vv_it ) = level;
						last_step++;
						vertex_count++;
			}}}
			if( last_step == 0 ) {
				if( continue_count < CONTINUE_COUNT_MAX )
					continue_count++;
				else break;
		}}
		return face_neighbors;
	}


/*
	-----------------------------------
	NEW VERTEX : COORDINATES AND NORMAL
	-----------------------------------
*/
	void compute_new_vertex(mesh_t& _m, typename mesh_t::FaceHandle face_handle) 
	{
		typedef std::set<typename mesh_t::VertexHandle> vhandles_t;
		typename vhandles_t::iterator v_it;
		typename mesh_t::FaceVertexIter fv_it = _m.fv_iter( face_handle );

		typename mesh_t::Point FNormal = _m.normal(face_handle);

		vhandles_t* face_neighborhood = face_neighbors( _m, face_handle );

		// prepare armadillo matrix, vectors
		arma::mat Gamma(10,10);
		Gamma.fill(0.0);

		arma::vec Omega(10);
		Omega.fill(0);

		// for each vertex in neighborhood
		for( v_it = face_neighborhood->begin(); v_it != face_neighborhood->end(); ++v_it )
		{
			typename mesh_t::Point
				Vcoords = _m.point( *v_it ),
				Ncoords = _m.normal( *v_it );

			double Vx =  Vcoords[0], Vy = Vcoords[1], Vz = Vcoords[2];
			double Nx =  Ncoords[0], Ny = Ncoords[1], Nz = Ncoords[2];

			double Vweight = VERTEX_WEIGHT_INIT;
			double Nweight = NORMAL_WEIGHT_INIT;

			for(int i=1; i < _m.property( vertex_level_, *v_it ); i++) {
				Vweight *= VERTEX_WEIGHT_FACTOR;
				Nweight *= NORMAL_WEIGHT_FACTOR;
			}

			arma::vec phi_i;
			phi_i << Vx*Vx << Vy*Vy << Vz*Vz << 2*Vx*Vy << 2*Vx*Vz << 2*Vy*Vz << 2*Vx << 2*Vy << 2*Vz << 1;

			arma::mat Psi_i;
			Psi_i	<< Vx*Vx 	<< 0 		<< 0 		<< Vx*Vy 		<< Vx*Vz 		<< 0 			<< Vx <<  0 <<  0 << 0 << arma::endr
					<< 0 		<< Vy*Vy 	<< 0 		<< Vx*Vy 		<< 0 			<< Vy*Vz 		<<  0 << Vy <<  0 << 0 << arma::endr
					<< 0 		<< 0 		<< Vz*Vz 	<< 0 			<< Vx*Vz 		<< Vy*Vz 		<<  0 <<  0 << Vz << 0 << arma::endr
					<< Vx*Vy 	<< Vx*Vy 	<< 0 		<< Vx*Vx+Vy*Vy 	<< Vy*Vz 		<< Vx*Vz 		<< Vy << Vx <<  0 << 0 << arma::endr
					<< Vx*Vz 	<< 0 		<< Vx*Vz 	<< Vy*Vz 		<< Vx*Vx+Vz*Vz 	<< Vx*Vy 		<< Vz <<  0 << Vx << 0 << arma::endr
					<< 0 		<< Vy*Vz	<< Vy*Vz 	<< Vx*Vz 		<< Vx*Vy 		<< Vy*Vy+Vz*Vz 	<<  0 << Vz << Vy << 0 << arma::endr
					<< Vx 		<< 0 		<< 0 		<< Vy 			<< Vz 			<< 0 			<<  1 <<  0 <<  0 << 0 << arma::endr
					<< 0 		<< Vy 		<< 0 		<< Vx 			<< 0 			<< Vz 			<<  0 <<  1 <<  0 << 0 << arma::endr
					<< 0 		<< 0 		<< Vz 		<< 0 			<< Vx 			<< Vy 			<<  0 <<  0 <<  1 << 0 << arma::endr
					<< 0 		<< 0 		<< 0 		<< 0 			<< 0 			<< 0 			<<  0 <<  0 <<  0 << 0 << arma::endr;

			Gamma += 2 * Vweight * phi_i * trans(phi_i) + 8 * Nweight * Psi_i;

			arma::vec Omega_i;
			Omega_i << Vx*Nx << Vy*Ny << Vz*Nz << Vx*Ny + Nx*Vy << Vx*Nz + Nx*Vz << Vy*Nz + Ny*Vz << Nx << Ny << Nz << 0;

			Omega += 4 * Nweight * Omega_i;
		}
		delete face_neighborhood;

		arma::vec solution = arma::solve(Gamma, Omega);
		_m.property( quadric_coeff_, face_handle ) = solution;


		// ---------------------------------------
		// PICKING THE NEW VERTEX FROM THE QUADRIC
		// ---------------------------------------

		//
		// Technique N°1 - projection of barycenter on quadric in the direction of triangle's normal
		//
		if( VERTEX_PICK_METHOD == 1 )
		{
			double a11 = solution(0), a22 = solution(1), a33 = solution(2);
			double a12 = solution(3), a13 = solution(4), a23 = solution(5);
			double a14 = solution(6), a24 = solution(7), a34 = solution(8);
			double a44 = solution(9);

			typename mesh_t::Point FCentroid = barycenterT( _m, face_handle );

			double FCx = FCentroid[0], FCy = FCentroid[1], FCz = FCentroid[2];
			double FNx = FNormal[0],   FNy = FNormal[1],   FNz = FNormal[2];

			// solve quadratic equation with coeff A, B, C
			double A;
				A  =    (a11 * FNx * FNx) + (a22 * FNy * FNy) + (a33 * FNz * FNz);
				A += 2*((a12 * FNx * FNy) + (a13 * FNx * FNz) + (a23 * FNy * FNz));

			double B;
				B  = (a11 * FCx * FNx) + (a22 * FCy * FNy) + (a33 * FCz * FNz);
				B +=  a12 * (FNx * FCy + FCx * FNy);
				B +=  a13 * (FNx * FCz + FCx * FNz);
				B +=  a23 * (FNy * FCz + FCy * FNz);
				B +=  a14 * FNx + a24 * FNy + a34 * FNz;
				B *= 2;

			double C;
				C  =    (a11 * FCx * FCx) + (a22 * FCy * FCy) + (a33 * FCz * FCz);
				C += 2* ((a12 * FCx * FCy) + (a13 * FCx * FCz) + (a23 * FCy * FCz));
				C += 2* (a14 * FCx + a24 * FCy + a34 * FCz) + a44;

			double discr = B*B - 4*A*C, t=0;
			
			if( discr < 0 ) t = 0;
			else 
			{
				double sqrt_discr = sqrt(discr);
				double t1 = (-B + sqrt_discr) / (2*A);
				double t2 = (-B - sqrt_discr) / (2*A);

				if( std::abs(t1) < std::abs(t2) )
					t = t1;
				else
					t = t2;
			}

			typename mesh_t::Point FNewVertex( FCx + t*FNx, FCy + t*FNy, FCz + t*FNz );
			_m.property( new_vertex_pos_, face_handle ) = FNewVertex;
			_m.property( new_vertex_normal_, face_handle ) = FNormal;
		} else {
		//
		// end Technique N°1
		//

		//
		// Technique N°2 - Foot Point
		//
			typename mesh_t::Point 
				bar_t = barycenterT( _m, face_handle ),
				new_vertex, new_normal;

			// compute new vertex
			new_vertex = footPoint(_m, face_handle);
			// choose new normal as n = gradient of f at new vertex + (new vertex - barycenter of T)
			new_normal = (new_vertex - bar_t) + quadric_gradient (_m, face_handle, new_vertex);

			// std::cout << new_normal << std::endl;

			if( scalar_product(new_normal, FNormal) < 0 )
				new_normal = -new_normal;

			// normalization?
			if( NORMALIZE_NORMALS )
				new_normal = new_normal.normalize();

			_m.property( new_vertex_pos_, face_handle ) = new_vertex;
			_m.property( new_vertex_normal_, face_handle ) = new_normal;
		}
		//
		// end Technique N°2
		//
	}




/*

·················································
··· FOOT POINT ALGORITHM FOR IMPLICIT SURFACE ···
·················································
*/

/*
	---------------------------
	FOOT POINT : MAIN PROCEDURE
	---------------------------
*/
	typename mesh_t::Point footPoint (mesh_t& _m, typename mesh_t::FaceHandle face_handle)
	{
		double error;
		typename mesh_t::Point p, p0, p1, q0, grad_p0, f1, f2, pp0;
		typename mesh_t::FaceVertexIter fv_it = _m.fv_iter( face_handle );

		p = barycenterT( _m, face_handle );
		p1 = surfacePoint( _m, face_handle, p );
		int step_count = 0;

		do {
			p0 = p1;
			pp0 = p - p0;

			grad_p0 = quadric_gradient( _m, face_handle, p0 );
			float coeff = scalar_product( pp0, grad_p0 ) / vector_length_sqr( grad_p0 );
			q0 = p - coeff * grad_p0;
			p1 = surfacePoint( _m, face_handle, q0);

			if( compute_error(q0, p0) > EPSILON )
			{
				f1 = q0 - p0;
				f2 = p1 - q0;

				float a0 = scalar_product( pp0, f1 );
				float a1 =  2 * scalar_product( pp0, f2 ) - vector_length_sqr( f1 );
				float a2 = -3 * scalar_product( f1, f2 );
				float a3 = -2 * vector_length_sqr( f2 );
				float alpha = 1 - (a0 + a1 + a2 + a3) / (a1 + 2*a2 + 3*a3);

				if( alpha > 0 && alpha < ALPHA_MAX )
				{
					q0 = p0 + alpha * f1 + (alpha * alpha) * f2;
					p1 = surfacePoint(_m, face_handle, q0);
			}}
			error = compute_error( p0, p1 );
		} while (error > EPSILON && ++step_count < 1000);
		return p1;
	}


/*
	--------------------------
	FOOT POINT : SURFACE POINT
	--------------------------
*/
	typename mesh_t::Point surfacePoint (mesh_t& _m, typename mesh_t::FaceHandle face_handle, typename mesh_t::Point q_start)
	{
		double error;
		typename mesh_t::Point q0, q1 = q_start, q0_grad;
		float q0_coeff;

		do {
			q0 = q1;
			q0_grad = quadric_gradient (_m, face_handle, q0);
			q0_coeff = quadric_value(_m, face_handle, q0) / vector_length_sqr(q0_grad);

			q1 = q0 - q0_coeff * q0_grad;

			error = compute_error( q0, q1 );
		} while (error > EPSILON);
		return q1;
	}


/*
	---------------------------------------------
	FOOT POINT : EVALUATION OF POLYNOMIAL f AT q0
	---------------------------------------------
*/
	float quadric_value (mesh_t& _m, typename mesh_t::FaceHandle face_handle, typename mesh_t::Point q0)
	{
		arma::vec sol = _m.property( quadric_coeff_, face_handle );

		double a11 = sol(0), a22 = sol(1), a33 = sol(2);
		double a12 = sol(3), a13 = sol(4), a23 = sol(5);
		double a14 = sol(6), a24 = sol(7), a34 = sol(8);
		double a44 = sol(9);

		double P_x = q0[0], P_y = q0[1], P_z = q0[2];

		return (float)
			a11*P_x*P_x +   a22*P_y*P_y +   a33*P_z*P_z +
		  2*a12*P_x*P_y + 2*a13*P_x*P_z + 2*a23*P_y*P_z +
		  2*a14*P_x     + 2*a24*P_y     + 2*a34*P_z     + a44;
	}


/*
	----------------------------------------------
	FOOT POINT : EVALUATION OF GRADIENT of f AT q0
	----------------------------------------------
*/
	typename mesh_t::Point quadric_gradient (mesh_t& _m, typename mesh_t::FaceHandle face_handle, typename mesh_t::Point q0)
	{
		arma::vec sol = _m.property( quadric_coeff_, face_handle );

		double a11 = sol(0), a22 = sol(1), a33 = sol(2);
		double a12 = sol(3), a13 = sol(4), a23 = sol(5);
		double a14 = sol(6), a24 = sol(7), a34 = sol(8);
		double a44 = sol(9);

		double P_x = q0[0], P_y = q0[1], P_z = q0[2];

		double grad_x = 2*a11*P_x + 2*a12*P_y + 2*a13*P_z + 2*a14;
		double grad_y = 2*a12*P_x + 2*a22*P_y + 2*a23*P_z + 2*a24;
		double grad_z = 2*a13*P_x + 2*a23*P_y + 2*a33*P_z + 2*a34;

		return typename mesh_t::Point(grad_x, grad_y, grad_z);
	}


/*
	---------------------------------
	FOOT POINT : COMPUTATION OF ERROR
	---------------------------------
*/
	double compute_error (typename mesh_t::Point& q0, typename mesh_t::Point& q1) 
	{
		typename mesh_t::Point q_dif = q0 - q1;
		return sqrt( vector_length_sqr( q_dif ) );
	}


/*
	-----------------------------
	FOOT POINT : LENGTH OF VECTOR
	-----------------------------
*/
	float vector_length_sqr (typename mesh_t::Point& Vec) 
	{
		double Vec_x = Vec[0];
		double Vec_y = Vec[1];
		double Vec_z = Vec[2];
		return (float) Vec_x*Vec_x + Vec_y*Vec_y + Vec_z*Vec_z;
	}

	float vector_length (typename mesh_t::Point& Vec) 
	{
		return sqrt( vector_length_sqr(Vec) );
	}

/*
	---------------------------
	FOOT POINT : SCALAR PRODUCT
	---------------------------
*/
	float scalar_product (typename mesh_t::Point& p, typename mesh_t::Point& q)
	{
		return (float) p[0]*q[0]+p[1]*q[1]+p[2]*q[2];
	}


/*
	---------------
	BARYCENTER of T
	---------------
 */
	typename mesh_t::Point barycenterT (mesh_t& _m, typename mesh_t::FaceHandle face_handle)
	{
		typename mesh_t::Point p;
		typename mesh_t::FaceVertexIter fv_it = _m.fv_iter( face_handle );
		p  = _m.point(*fv_it);
		p += _m.point(*(++fv_it));
		p += _m.point(*(++fv_it));
		p *= 1.00 / 3.00;
		return p;
	}


/*
	-------------------------------
	TRIANGLE AREA - Heron's formula
	-------------------------------
 */
	float triangleArea (mesh_t& _m, typename mesh_t::FaceHandle face_handle)
	{
		typename mesh_t::Point A, B, C;
		float a, b, c, p;
		typename mesh_t::FaceVertexIter fv_it = _m.fv_iter( face_handle );

		A = _m.point(*fv_it);
		B = _m.point(*(++fv_it));
		C = _m.point(*(++fv_it));

		a = vector_length(B - C);
		b = vector_length(A - C);
		c = vector_length(A - B);
		p = (a + b + c) * 0.5;
		
		return (float) sqrt( p * (p-a) * (p-b) * (p-c) );
	}
};


//=============================================================================
} // END_NS_UNIFORM
} // END_NS_SUBDIVIDER
} // END_NS_OPENMESH
//=============================================================================
#endif // OPENMESH_SUBDIVIDER_UNIFORM_QUADRATICFITTINGT_HH defined
//=============================================================================
