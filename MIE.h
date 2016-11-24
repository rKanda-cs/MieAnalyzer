//=============================================================================                              
// Copyright (C) 2013 Nyam-Erdene Odontsengel.
// Last updated by on Nov, 2013.
//=============================================================================

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <complex>
#include <iostream>
using namespace std;

complex<double> I = complex<double>(0, 1);

// ### choise
//const int ITERATION = 100;
const int ITERATION =  20;

// bessel function of first kind
complex<double> Bessel1(int order, double x){
	switch( order ){
		case 0:
			return _j0( x );
		case 1:
			return _j1( x );
		default:
			return _jn( order, x );
	}
}

// bessel function of second kind
complex<double> Bessel2(int order, double x){
	switch( order ){
		case 0:
			return _y0( x );
		case 1:
			return _y1( x );
		default:
			return _yn( order, x );
	}
}

// hankel function of first kind
complex<double> Hankel1(int order, double x){
	return Bessel1(order, x) + I * Bessel2(order, x);
}

// hankel function of second kind
complex<double> Hankel2(int order, double x){
	return Bessel1(order, x) - I * Bessel2(order, x);
}

// differential of Bessel function of first kind
complex<double> d_Bessel1(int order, double x){
	return Bessel1(order-1, x) - Bessel1(order, x).operator *= (order / x);
}
// differential of Hankel function of first kind
complex<double> d_Hankel1(int order, double x){
	return Hankel1(order-1, x) - Hankel1(order, x).operator *= (order / x);
}
// coefficient an for Mie scattering
complex<double> a_order(int order, double x, double index){
	return ( d_Bessel1(order, index*x) * Bessel1(order, x) - index * Bessel1(order, index*x) * d_Bessel1(order, x) )
		/  ( d_Bessel1(order, index*x) * Hankel1(order, x) - index * Bessel1(order, index*x) * d_Hankel1(order, x) );
}

// coefficient cn for Mie scattering
complex<double> c_order(int order, double x, double index){
	return ( Bessel1(order, x) - a_order(order, x, index) * Hankel1(order, x) )
		/  Bessel1(order, index*x);
}

// calculate a_order by iteration
/*
complex<double> A_order(int order, double index, double x){
	double order_c = x + 4.05*pow(x, 1/3) + 2;
	int order_index_x = max(order_c, abs(index*x)) + 15;
	cout << "order_c=" << order_c << " index*x=" << index*x <<" order_index_x=" << order_index_x << endl;

	complex<double> A_order_index_x = complex<double>(0, 0);

	for(double i=order_index_x-1; i>=order; i=i-1){
		A_order_index_x = i/(index*x) - complex<double>(1, 0)/(A_order_index_x+(i+1)/(index*x));
		cout << "A_order_index_x[" << i << "]=" << A_order_index_x << endl; 
	}

	return A_order_index_x;
}
*/

////////////////////////////////////////////////////////////////////////////////////
// (r,th: polar coordinate, k: wavenumber, radius: radius, index: refractive index) 
// polar coordinate
// incident field
complex<double> Er_0(double k, double r, double th){
	return exp(I * k * r * cos(th)) * sin(th);
}
complex<double> Eth_0(double k, double r, double th){
	return exp(I * k * r * cos(th)) * cos(th);
}

// internal scattered field
complex<double> Er_i(double index, double radius, double k, double r, double th){
	complex<double> E = 0;
	int order;

	for(order=1; order<=ITERATION; order++){
		E += pow(I, order) * c_order(order, k * radius, index) * complex<double>(order, 0) * Bessel1(order, index * k * r) * sin(order * th);
	}
	E = - I * complex<double>(2/(index*index*k*r), 0) * E - Er_0(k, r, th);

	return E;
}
complex<double> Eth_i(double index, double radius, double k, double r, double th){
	complex<double> E = 0;
	int order;

	for(order=1; order<=ITERATION; order++){
		E += pow(I, order) * c_order(order, k * radius, index) * d_Bessel1(order, index * k * r) * cos(order * th);
	}
	E = - I * complex<double>(1/index, 0) * c_order(0, k*radius, index) * d_Bessel1(0, index * k * r) - I * complex<double>(2/index, 0) * E - Eth_0(k, r, th);

	return E;
}

// external scattered field
complex<double> Er_s(double index, double radius, double k, double r, double th){
	complex<double> E = 0;
	int order;

	for(order=1; order<=ITERATION; order++){
		E += pow(I, order) * a_order(order, k * radius, index) * complex<double>(order, 0) * Hankel1(order, k * r) * sin(order * th);
	}
	E = I * complex<double>(2/(k*r), 0) * E;

	return E;
}
complex<double> Eth_s(double index, double radius, double k, double r, double th){
	complex<double> E = 0;
	int order;

	for(order=1; order<=ITERATION; order++){
		E += pow(I, order) * a_order(order, k * radius, index) * d_Hankel1(order, k * r) * cos(order * th);
	}
	E = I * ( a_order(0, k*radius, index) * d_Hankel1(0, k * r) + complex<double>(2, 0) * E );

	return E;
}

////////////////////////////////////////////////////////////////////////////////////
// Cartesian coordinate
// (r,th: polar coordinate, k: wavenumber, radius: radius, index: refractive index)
// internal scattered field
complex<double> TEx_scattered_internal( double r, double th, double k, double radius, double index ){
	return Er_i(index, radius, k, r, th) * cos(th) - Eth_i(index, radius, k, r, th) * sin(th);
}
complex<double> TEy_scattered_internal( double r, double th, double k, double radius, double index ){
	return Er_i(index, radius, k, r, th) * sin(th) + Eth_i(index, radius, k, r, th) * cos(th);
}

// external scattered field
complex<double> TEx_scattered_external( double r, double th, double k, double radius, double index ){
	return Er_s(index, radius, k, r, th) * cos(th) - Eth_s(index, radius, k, r, th) * sin(th);
}
complex<double> TEy_scattered_external( double r, double th, double k, double radius, double index ){
	return Er_s(index, radius, k, r, th) * sin(th) + Eth_s(index, radius, k, r, th) * cos(th);
}
