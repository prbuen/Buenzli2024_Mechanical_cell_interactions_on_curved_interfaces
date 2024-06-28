/*
 * Author: Pascal R Buenzli, 2022-2024
 */

module model_inputs;

import std.stdio;
import std.exception;
import std.math : sqrt,cos,sin,atan2,PI,lround;
import std.numeric : dotProduct, euclideanDistance, findRoot;
import std.algorithm: max;

import utils;
import main;

alias Real=double;
alias NodeArrayT=Node[];

// ===== choose model features, versions =====
// version = Model0; // circular interface
// version = Model2; // cross
// version = Model2b; // cross with 7 cells evenly distributed except for one boundary
// version = Model2c; // cross with 8 cells and deleting one node => 7 cells
version = Model4; // open sin curve with stressed or stressfree steady state

version = hookean_force; // f(l)=k(l-a)
// version = lineardiff_force; // f(l)=a^2 k(1/a-1/l), so that f'(a)=k in both models
// version = porousdiff_force; // f(l) = (a^3 k/2)(1/a^2 - 1/l^2) so that f'(a)=k and D(q)=D0 q, with D0 = a^3 k/eta

version = straight_springs;
// version = curved_springs;

version = rescale_spring_params; // with number of nodes, for continuum limit

// version = closed_curve; // automatically set in each model
// version = open_curve;

class ModelInputs_base {
public:
	string model = "";

	// export compilation switches
	mixin(DefVersion!"hookean_force");
	mixin(DefVersion!"lineardiff_force");
	mixin(DefVersion!"porousdiff_force");
	mixin(DefVersion!"straight_springs");
	mixin(DefVersion!"curved_springs");
	mixin(DefVersion!"rescale_spring_params");
	mixin(DefVersion!"integral_simpson");
	mixin(DefVersion!"closed_curve");
	mixin(DefVersion!"open_curve");

	// data storage
	Data data;

	// time parameters
	Real t_end;
	Real dt;
	ulong output_every_nth_step;
	Real dt_frame; // = dt*output_every_nth_step (convenience)

	// interface parameters
	Real R; // circular curve radius
	Real u_max; // interface parameter max
	Real s_max; // length of interface (arclength max)

	// cell parameters
	Real k; // spring constant (scaling factor)
	Real a; // spring resting length
	Real eta; // drag coefficient

	size_t N; // number of cells
	size_t m; // number of springs per cell
	size_t M; // number of springs = N*m (convenience)
	size_t M_nodes; // number of spring boundaries: M (loop) or M+1 (open curve)

	/* force law scalar in the tangential direction. Depends on signed distance l.
	   Positive in extension, negative in compression
	   Independent of straight vs curved springs
	 */
	Real force_tg(Real l)
	{
		static if(hookean_force)
			return k*(l-a);

		else static if(lineardiff_force)
			return a^^2*k*(1./a - 1./l);

		else static if(porousdiff_force)
			return 0.5*a^^3*k*(1./a^^2 - 1./l^^2);

		enforce(0, "Not implemented");
		return 0;
	}

	Real[2] r(Real u) { return [0,0]; }
	Real[2] tau(Real u) { return [0,0]; }
	Real arclength(Real u) { return 0; }
	Real kappa(Real u) { return 0; }
	Real u(Real s)
	{
		auto root = findRoot( (Real u) => arclength(u)-s, 0.0, u_max, -s, s_max-s);
		return root[0];
	}


	// For correct calculations of arclength distance between nodes, node indices must be ordered by increasing arclength parameter
	void initial_springs_props(ref NodeArrayT nodes)
	{
		enforce(nodes.length==M_nodes && M==m*N && M != 0 && M_nodes != 0, "No nodes");

		// Set springs according to arclength first then calculate the corresponding u.
		
		// set cell boundaries: 0, m, 2m, 3m, ...
		auto offset = 1e-9*s_max;
		static if(open_curve)
		{
			offset = 0;
			nodes[N*m].s = s_max;
		}
		foreach(i; 0 .. N)
			nodes[i*m].s = offset + i*s_max/N; // overall offset so PBC boundary is within last spring (otherwise, there are issues with connecting last and first)

		// offset one of the cell boundaries
		size_t i0 = N/2;// N/2; 1; N-1;
		nodes[i0*m].s += 0.5*s_max/N;

		// i0 = N/4;
		// nodes[i0*m].s += 0.4*s_max/N;
		// nodes[(i0+1)*m].s += 0.4*s_max/N;
		
		// divide each cell into m equidistant springs: 1,2,..,m-1, m+1,m+2,...
		foreach(i; 0 .. N)
		{
			size_t ip1 = (i==N-1? 0 : i+1);
			foreach(j; 1 .. m)
				nodes[i*m+j].s = nodes[i*m].s + j*s_PBC((nodes[ip1*m].s - nodes[i*m].s))/m;
		}

		// calculate u, pos vectors r, and (true) tangents from the position along circle
		foreach(ref n; nodes)
		{
			n.u = u(n.s);
			n.r = r(n.u);
			n.tau = tau(n.u);
		}

		// calculate length (1/density) and stress (between current and previous)
		foreach(i, ref n; nodes)
		{
			size_t im1 = (i==0? nodes.length-1 : i-1);
			static if(straight_springs)
				n.length = euclideanDistance(n.r[], nodes[im1].r[]);
			else static if(curved_springs)
				n.length = s_PBC(n.s - nodes[im1].s);
			else
				enforce(0, "not implemented");
			n.stress = force_tg(n.length)/(k*a); // F/(k a) = sigma/E, E=Young's modulus WARNING: normalised by resting length a, not current length!
		}
	}
	
	// parameter wraparound. Suitable for u and for u2-u1
	// Not safe if |u| > 2*u_max: best to use remainder(u, u_max) = sign(u/u_max)*(u-n*u_max) or remquo(u, u_max)
	pragma(inline): Real u_PBC(Real u)
	{
		return u<0? u+u_max : (u >= u_max? u-u_max : u);
	}
	
	// arclength wraparound. Suitable for s and for s2-s1
	// Not safe if |ss| > 2*s_max: best to use remainder(s, s_max) = sign(s/s_max)*(s-n*s_max) or remquo(s, s_max)
	pragma(inline): Real s_PBC(Real s)
	{
		return s<0? s+s_max : (s >= s_max? s-s_max : s);
	}
}

version(Model0) {
	version = closed_curve;
class ModelInputs : ModelInputs_base {
public:

	this()
	{
		model = "Model0";
		data = Data("datadir/springs_0.dat");
		data.erase_datadir();
	
		// time parameters
		t_end = 2;
		dt = 0.001;
		output_every_nth_step = 1;
		dt_frame = dt*output_every_nth_step;

		// interface parameters
		R = 1;
		u_max = 2*PI;
		s_max = arclength(u_max);

		// cell parameters (may be rescaled by m)
		k = 1;
		a = 1;
		eta = 1;
	
		N = 8;
		m = 4;
		M = m*N;
		M_nodes = M; static if (ModelInputs.open_curve) { M_nodes++; }

		static if(curved_springs)
		{
			// Uniformly distributed springs have no stress gradient
			a = s_max/M;
			// a = 0; // always in extension
		}
		else static if(straight_springs)
		{
			// For a general interface, we can't initiate positions without stress gradient with straight springs. What we can do is initiate with uniform arclength then let the system evolve to equilibrium.

			// For now, we simply initiate with uniform arclength
			a = s_max/M;
			
			// a = 0; // always in extension

			// For M springs on a circle the steady state is a regular M-gon (see https://en.wikipedia.org/wiki/Regular_polygon, circumradius)
			a = 2*R*sin(PI/M); // NB: for M->inf this gives 2*pi*R/M, ok.
		}
		else
			enforce(0, "not implemented");

		static if(rescale_spring_params)
		{
			// rescale k and eta (a is already rescaled above)
			k *= m;
			eta /= m;
		}
		
	}
	

	// Interface parametrisation and corresponding (true) tangent
	// u is a general parameter, s is arclength. 
	override Real[2] r(Real u)
	{
		return [R*cos(u), R*sin(u)];
	}
	override Real[2] tau(Real u)
	{
		return [-sin(u), cos(u)];
	}
	override Real arclength(Real u)
	{
		return u*R;
	}
	override Real u(Real s)
	{
		return s/R; 
	}
	override Real kappa(Real u)
	{
		return -R;
	}
}
}

version(Model2) {
	version = closed_curve;
class ModelInputs : ModelInputs_base {
public:
	this()
	{
		model = "Model2";
		data = Data("datadir/springs_0.dat");
		data.erase_datadir();

		// time parameters
		// t_end = 100.;
		// output_every_nth_step = 100;
		t_end = 4.;
		output_every_nth_step = 25;
		dt = 0.001;
		dt_frame = dt*output_every_nth_step;

		// interface parameters
		R = 1.;
		u_max = 2*PI;
		s_max = arclength(u_max);

		// cell parameters (may be rescaled)
		N = 8;
		m = 1;
		M = N*m;
		M_nodes = M; static if (ModelInputs.open_curve) { M_nodes++; }

		k = 1;
		eta = 1;
		static if(rescale_spring_params)
		{
			// rescale k and eta (a is already rescaled above)
			k *= m;
			eta /= m;
		}
			
	
		static if(curved_springs)
		{
			// Uniformly distributed springs have no stress gradient
			a = s_max/M;
			// a = 0; // always in extension
		}
		else static if(straight_springs)
		{
			// For a general interface, we can't initiate positions without stress gradient with straight springs. What we can do is initiate with uniform arclength then let the system evolve to equilibrium.

			// For now, we simply initiate with uniform arclength
			if (m==1)
				a = 0.736813; // empirical with N=8,m=1,t=20 (no change with t=100)
			else if (m==4)
				a= 0.196175333333333; // empirical with N=8,m=4,t=20
			else if (m==8)
				a = 0.0994770769230769; // empirical with N=8,m=8,t=30, rel stdev(a)=0.0229
			else
				a = s_max/M;

			// a = s_max/M; // used for determination of a for stree-free equilibrium
		}
		else
			enforce(0, "not implemented");
	}


	// Interface parametrisation and corresponding (true) tangent
	// u is a general parameter, s is arclength. 
	override Real[2] r(Real u)
	{
		Real r_theta = R*( (cos(u))^^4 + (sin(u))^^4 );
		return [r_theta*cos(u), r_theta*sin(u)];
	}
	override Real[2] tau(Real u)
	{
		Real r_theta = R*( (cos(u))^^4 + (sin(u))^^4 );
		Real rp_theta = 4*R*( -sin(u)*(cos(u))^^3 + cos(u)*(sin(u))^^3 );		
		Real[2] rprime = rp_theta*[cos(u), sin(u)] + r_theta*[-sin(u), cos(u)];
		rprime[] /= norm(rprime[]);
		return rprime;
	}
	override Real arclength(Real u)
	{
		ulong n = 1001;
		Real[] v = new Real[n];
		foreach(i, ref vi; v)
		{
			auto ui = i*u/n;
			Real r_theta = R*( (cos(ui))^^4 + (sin(ui))^^4 );
			Real rp_theta = 4*R*( -sin(ui)*(cos(ui))^^3 + cos(ui)*(sin(ui))^^3 );		
			Real[2] rprime = rp_theta*[cos(ui), sin(ui)] + r_theta*[-sin(ui), cos(ui)];
			vi = norm(rprime[]);
			// stderr.writeln("v_i = ", vi);
		}
		return simpson(u-0, v);
	}	
}
}


version(Model2b) {
	version = closed_curve;
class ModelInputs : ModelInputs_base {
public:
	this()
	{
		model = "Model2b";
		data = Data("datadir/springs_0.dat");
		data.erase_datadir();

		// time parameters
		// t_end = 100.;
		// output_every_nth_step = 100;
		t_end = 4.;
		output_every_nth_step = 25;
		dt = 0.001;
		dt_frame = dt*output_every_nth_step;

		// interface parameters
		R = 1.;
		u_max = 2*PI;
		s_max = arclength(u_max);

		// cell parameters (may be rescaled)
		N = 7;
		m = 8;
		M = N*m;
		M_nodes = M; static if (ModelInputs.open_curve) { M_nodes++; }

		k = 1;
		eta = 1;
		static if(rescale_spring_params)
		{
			// rescale k and eta (a is already rescaled above)
			k *= m;
			eta /= m;
		}
			
	
		static if(curved_springs)
		{
			// Uniformly distributed springs have no stress gradient
			a = s_max/M;
			// a = 0; // always in extension
		}
		else static if(straight_springs)
		{
			// For a general interface, we can't initiate positions without stress gradient with straight springs. What we can do is initiate with uniform arclength then let the system evolve to equilibrium.

			// For now, we simply initiate with uniform arclength
			if (m==1)
				// a = 0.736813; // empirical with N=8,m=1,t=20 (no change with t=100)
				a = 0.7364893; // empirical with N=7,m=1,t=20 (no change with t=100)
			else if (m==4)
				a= 0.196175333333333; // empirical with N=8,m=4,t=20
			else if (m==8)
				// a = 0.0994770769230769; // empirical with N=8,m=8,t=30, rel stdev(a)=0.0229
				a = 0.113400859649123; // empirical with N=7,m=8,t=100, rel 
			else
				a = s_max/M;

			// a = s_max/M; // used for determination of a for stree-free equilibrium
		}
		else
			enforce(0, "not implemented");
	}


	// Interface parametrisation and corresponding (true) tangent
	// u is a general parameter, s is arclength. 
	override Real[2] r(Real u)
	{
		Real r_theta = R*( (cos(u))^^4 + (sin(u))^^4 );
		return [r_theta*cos(u), r_theta*sin(u)];
	}
	override Real[2] tau(Real u)
	{
		Real r_theta = R*( (cos(u))^^4 + (sin(u))^^4 );
		Real rp_theta = 4*R*( -sin(u)*(cos(u))^^3 + cos(u)*(sin(u))^^3 );		
		Real[2] rprime = rp_theta*[cos(u), sin(u)] + r_theta*[-sin(u), cos(u)];
		rprime[] /= norm(rprime[]);
		return rprime;
	}
	override Real arclength(Real u)
	{
		ulong n = 1001;
		Real[] v = new Real[n];
		foreach(i, ref vi; v)
		{
			auto ui = i*u/n;
			Real r_theta = R*( (cos(ui))^^4 + (sin(ui))^^4 );
			Real rp_theta = 4*R*( -sin(ui)*(cos(ui))^^3 + cos(ui)*(sin(ui))^^3 );		
			Real[2] rprime = rp_theta*[cos(ui), sin(ui)] + r_theta*[-sin(ui), cos(ui)];
			vi = norm(rprime[]);
			// stderr.writeln("v_i = ", vi);
		}
		return simpson(u-0, v);
	}


		// For correct calculations of arclength distance between nodes, node indices must be ordered by increasing arclength parameter
	// Here we want to remove the 3rd node out of 8, so we have 7 left.
	override void initial_springs_props(ref NodeArrayT nodes)
	{
		enforce(nodes.length==M_nodes && M==m*N && M != 0 && M_nodes != 0, "No nodes");

		// Set springs according to arclength first then calculate the corresponding u.
		
		// set cell boundaries: 0, m, 2m, 3m, ...
		auto offset = 1e-9*s_max;
		static if(open_curve)
		{
			offset = 0;
			nodes[N*m].s = s_max;
		}
		foreach(i; 0 .. N)
			nodes[i*m].s = offset + i*s_max/N; // overall offset so PBC boundary is within last spring (otherwise, there are issues with connecting last and first)

		// offset one of the cell boundaries
		size_t i0 = N/2;// N/2; 1; N-1;
		nodes[i0*m].s += 0.5*s_max/N;

		// i0 = N/4;
		// nodes[i0*m].s += 0.4*s_max/N;
		// nodes[(i0+1)*m].s += 0.4*s_max/N;
		
		// divide each cell into m equidistant springs: 1,2,..,m-1, m+1,m+2,...
		foreach(i; 0 .. N)
		{
			size_t ip1 = (i==N-1? 0 : i+1);
			foreach(j; 1 .. m)
				nodes[i*m+j].s = nodes[i*m].s + j*s_PBC((nodes[ip1*m].s - nodes[i*m].s))/m;
		}

		// calculate u, pos vectors r, and (true) tangents from the position along circle
		foreach(ref n; nodes)
		{
			n.u = u(n.s);
			n.r = r(n.u);
			n.tau = tau(n.u);
		}

		// calculate length (1/density) and stress (between current and previous)
		foreach(i, ref n; nodes)
		{
			size_t im1 = (i==0? nodes.length-1 : i-1);
			static if(straight_springs)
				n.length = euclideanDistance(n.r[], nodes[im1].r[]);
			else static if(curved_springs)
				n.length = s_PBC(n.s - nodes[im1].s);
			else
				enforce(0, "not implemented");
			n.stress = force_tg(n.length)/(k*a); // F/(k a) = sigma/E, E=Young's modulus WARNING: normalised by resting length a, not current length!
		}
	}

}
}


version(Model2c) {
	version = closed_curve;
class ModelInputs : ModelInputs_base {
public:
	this()
	{
		model = "Model2c";
		data = Data("datadir/springs_0.dat");
		data.erase_datadir();

		// time parameters
		// t_end = 100.;
		// output_every_nth_step = 100;
		t_end = 4.;
		output_every_nth_step = 25;
		dt = 0.001;
		dt_frame = dt*output_every_nth_step;

		// interface parameters
		R = 1.;
		u_max = 2*PI;
		s_max = arclength(u_max);

		// cell parameters (may be rescaled)
		N = 7;
		m = 8;
		M = N*m;
		M_nodes = M; static if (ModelInputs.open_curve) { M_nodes++; }

		k = 1;
		eta = 1;
		static if(rescale_spring_params)
		{
			// rescale k and eta (a is already rescaled above)
			k *= m;
			eta /= m;
		}
			
	
		static if(curved_springs)
		{
			// Uniformly distributed springs have no stress gradient
			a = s_max/M;
			// a = 0; // always in extension
		}
		else static if(straight_springs)
		{
			// For a general interface, we can't initiate positions without stress gradient with straight springs. What we can do is initiate with uniform arclength then let the system evolve to equilibrium.

			// For now, we simply initiate with uniform arclength
			if (m==1)
				// a = 0.736813; // empirical with N=8,m=1,t=20 (no change with t=100)
				a = 0.7364893; // empirical with N=7,m=1,t=20 (no change with t=100)
			else if (m==4)
				a= 0.196175333333333; // empirical with N=8,m=4,t=20
			else if (m==8)
				// a = 0.0994770769230769; // empirical with N=8,m=8,t=30, rel stdev(a)=0.0229
				a = 0.113400859649123; // empirical with N=7,m=8,t=100, rel 
			else
				a = s_max/M;

			// a = s_max/M; // used for determination of a for stree-free equilibrium
		}
		else
			enforce(0, "not implemented");
	}


	// Interface parametrisation and corresponding (true) tangent
	// u is a general parameter, s is arclength. 
	override Real[2] r(Real u)
	{
		Real r_theta = R*( (cos(u))^^4 + (sin(u))^^4 );
		return [r_theta*cos(u), r_theta*sin(u)];
	}
	override Real[2] tau(Real u)
	{
		Real r_theta = R*( (cos(u))^^4 + (sin(u))^^4 );
		Real rp_theta = 4*R*( -sin(u)*(cos(u))^^3 + cos(u)*(sin(u))^^3 );		
		Real[2] rprime = rp_theta*[cos(u), sin(u)] + r_theta*[-sin(u), cos(u)];
		rprime[] /= norm(rprime[]);
		return rprime;
	}
	override Real arclength(Real u)
	{
		ulong n = 1001;
		Real[] v = new Real[n];
		foreach(i, ref vi; v)
		{
			auto ui = i*u/n;
			Real r_theta = R*( (cos(ui))^^4 + (sin(ui))^^4 );
			Real rp_theta = 4*R*( -sin(ui)*(cos(ui))^^3 + cos(ui)*(sin(ui))^^3 );		
			Real[2] rprime = rp_theta*[cos(ui), sin(ui)] + r_theta*[-sin(ui), cos(ui)];
			vi = norm(rprime[]);
			// stderr.writeln("v_i = ", vi);
		}
		return simpson(u-0, v);
	}


		// For correct calculations of arclength distance between nodes, node indices must be ordered by increasing arclength parameter
	// Here we want to remove the 3rd node out of 8, so we have 7 left.
	override void initial_springs_props(ref NodeArrayT nodes)
	{
		enforce(nodes.length==M_nodes && M==m*N && M != 0 && M_nodes != 0, "No nodes");

		// Set springs according to arclength first then calculate the corresponding u.

		// Positions are based on 8 springs, and removing one node
		
		// set cell boundaries: 0, m, 2m, 3m, ...
		auto offset = 1e-9*s_max;
		static if(open_curve)
		{
			offset = 0;
			nodes[(N+1)*m].s = s_max;
		}
		foreach(i; 0 .. N+1)
		{
			if (i<3)
				nodes[i*m].s = offset + i*s_max/(N+1); // overall offset so PBC boundary is within last spring (otherwise, there are issues with connecting last and first)
			else
				nodes[(i-1)*m].s = offset + i*s_max/(N+1);
		}
		// offset one of the cell boundaries
		size_t i0 = (N+1)/2-1;// N/2; 1; N-1;
		nodes[i0*m].s += 0.5*s_max/(N+1);

		// i0 = N/4;
		// nodes[i0*m].s += 0.4*s_max/N;
		// nodes[(i0+1)*m].s += 0.4*s_max/N;
		
		// divide each cell into m equidistant springs: 1,2,..,m-1, m+1,m+2,...
		foreach(i; 0 .. N)
		{
			size_t ip1 = (i==N-1? 0 : i+1);
			foreach(j; 1 .. m)
				nodes[i*m+j].s = nodes[i*m].s + j*s_PBC((nodes[ip1*m].s - nodes[i*m].s))/m;
		}

		// calculate u, pos vectors r, and (true) tangents from the position along circle
		foreach(ref n; nodes)
		{
			n.u = u(n.s);
			n.r = r(n.u);
			n.tau = tau(n.u);
		}

		// calculate length (1/density) and stress (between current and previous)
		foreach(i, ref n; nodes)
		{
			size_t im1 = (i==0? nodes.length-1 : i-1);
			static if(straight_springs)
				n.length = euclideanDistance(n.r[], nodes[im1].r[]);
			else static if(curved_springs)
				n.length = s_PBC(n.s - nodes[im1].s);
			else
				enforce(0, "not implemented");
			n.stress = force_tg(n.length)/(k*a); // F/(k a) = sigma/E, E=Young's modulus WARNING: normalised by resting length a, not current length!
		}
	}

}
}

version(Model4) {
	version = open_curve;
class ModelInputs : ModelInputs_base {
public:
	this()
	{
		model = "Model4";
		data = Data("datadir/springs_0.dat");
		data.erase_datadir();
		
		// time parameters
		t_end = 4;
		dt = 0.001;
		output_every_nth_step = 10;
		dt_frame = dt*output_every_nth_step;

		// interface parameters
		R = 0.8; // amplitude
		u_max = 2*PI;
		s_max = arclength(u_max);

		// cell parameters (may be rescaled)
		N = 4;
		m = 4;
		M = N*m;
		M_nodes = M; static if (ModelInputs.open_curve) { M_nodes++; }

		k = 1;
		eta = 1;
		static if(rescale_spring_params)
		{
			// rescale k and eta (a is already rescaled above)
			k *= m;
			eta /= m;
		}
	
		static if(curved_springs)
		{
			// Uniformly distributed springs have no stress gradient
			a = s_max/M;
			// a = 0; // always in extension
		}
		else static if(straight_springs)
		{
			// For a general interface, we can't initiate positions without stress gradient with straight springs. What we can do is initiate with uniform arclength then let the system evolve to equilibrium.

			// For now, we simply initiate with uniform arclength
			a = s_max/M;
			
			// a = 0; // always in extension
		}
		else
			enforce(0, "not implemented");

		a *= 1; // stressfree steady state
		// a *= 0.5; // stressed steady state in extension
		// a *= 2; // stressed steady state in compression
	}


	// Interface parametrisation and corresponding (true) tangent
	// u is a general parameter, s is arclength. 
	override Real[2] r(Real u)
	{
		return [u, R*sin(u)];
	}
	override Real[2] tau(Real u)
	{
		Real[2] rprime = [1, R*cos(u)];
		rprime[] /= norm(rprime[]);
		return rprime;
	}
	override Real arclength(Real u)
	{
		size_t n = 1001;
		Real[] v = new Real[n];
		foreach(i, ref vi; v)
		{
			Real ui = i*u/n;
			Real[2] rprime = [1, R*cos(ui)];
			vi = norm(rprime[]);
			// stderr.writeln("v_i = ", vi);
		}
		return simpson(u-0, v);
	}

	// curvature: r(u)=(x(u), y(u)) then kappa(u) = (y_uu x_u - x_uu y_u)/(x_u^2+y_u^2)^(3/2)
	override Real kappa(Real u)
	{
		auto x_u = 1;
		auto x_uu = 0;
		auto y_u = R*cos(u);
		auto y_uu = -R*sin(u);

		return (y_uu*x_u - x_uu*y_u)/(x_u^^2 + y_u^^2)^^(3/2);
	}
}
}
