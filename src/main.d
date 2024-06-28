/*
 * Author: Pascal R Buenzli, 2022-2024
 *
 * Description:
 * Buenzli PR, Kuba S, Murphy RJ, Simpson MJ (2024) Mechanical cell interactions on curved interfaces. Preprint available at https://arxiv.org/TBA
 *
 * Requirements:
 * D compiler: ldc2 or dmd, see https://wiki.dlang.org/LDC
 * For visualisation: Gnuplot > 5.4.x and Python
 *
 * Compilation:
 * Edit model_inputs.d to select model and parameters sets, then compile:
 * 	ldc2 src/*.d -od=build
 * or (optimised):
 * 	ldc2 src/*.d -od=build -O3 -release -mcpu=native -ffast-math
 * 
 * Execution:
 * ./main
 * => data created in `datadir/*.dat` (filenames are suffixed with time frame number)
 * 
 * Visualisation:
 * Edit snapshot.py, trajectories.py, or profile.py (e.g. to set frame=70) and pipe into gnuplot (alternatively, save python output to file and load in gnuplot)
 *   python snapshot.py | gnuplot
 * => snapshot.pdf
 * or
 *   python trajectories.py | gnuplot
 * => trajectories.pdf
 * or
 *   python profiles.py | gnuplot
 * => profiles.pdf
 */

module main;

import std.stdio;
import std.exception;
import std.range: iota;
import std.math: sqrt,cos,sin,atan2,PI,lround;
import std.numeric: dotProduct, euclideanDistance;

import utils;
import model_inputs: Real, ModelInputs;

struct Node {
	Real[2] r; // position
	Real[2] tau; // unit tangent (positively oriented, Jordan curve)
	Real u; // interface parameter
	Real s; // corresponding arc length parameter
	Real length; // length between current and _previous_ node
	Real stress; // F/(k*a)=sigma/E (E: Young's modulus), between current and _previous_ node, a = resting length (not current length!)
}

class Model {
protected:
	// model inputs
	ModelInputsWrap!(ModelInputs) _p;
	alias _p this; // make ModelInputs members part of Model

	// state variables
	Real t;
	ulong iter; // step() iteration
	ulong frame_iter; // write_state() iteration	
	Node[] nodes;
	Node[] old_nodes; // for time stepping

public:
	this() // initialise
	{
		_p = new ModelInputsWrap!(ModelInputs)(); // implicitly defines p

		// write parameters to stderr, and keep a copy in datadir
		stderr.writeln(_p.params() ~ "\n");
		stderr.writeln("#### Data written to ", data.datadir, " ####");
		stderr.writeln("#### Call replay on ", data.filepath(data.name_pattern, data.init_frame_id), " ####");

		string metadata = _p.params("# ", ["description", "model_version", "firstdatafile", "data", "lambda"] ~ Functions); // comment out strings and functions, and parameter lambda (=python keyword)
		auto metadata_fp = data.open("metadata.dat");
		metadata_fp.write(metadata);
		metadata_fp.flush();

		// state variables
		t = 0;
		iter = 0;

		nodes = new Node[M_nodes];
		old_nodes = new Node[M_nodes];
		
		initial_springs_props(nodes);
	}


	void write_state()
	{
		auto springs_data = data.open("springs", frame_iter);
		springs_data.writeln("# u, arclength, x, y, tau_1, tau_2, length, stress, curvature");

		// spring and node data
		// NB: for series of data like x_i, y_i, f_i, gnuplot colours the segment i->i+1 with f_i+1, so we prepare f_i so it corresponds to the length or force along i-1 -> i.
		foreach(i, n; nodes)
			with(n)
			{
				springs_data.writeln(u, "\t", s, "\t", r[0], "\t", r[1], "\t", tau[0], "\t", tau[1], "\t", length, "\t", stress, "\t", kappa(u));
			}

		// wraparound for spring data and ease of plotting
		static if(closed_curve)
		{
			with(nodes[0])
			{
				springs_data.writeln(u, "\t", s, "\t", r[0], "\t", r[1], "\t", tau[0], "\t", tau[1], "\t", length, "\t", stress, "\t", kappa(u));			
			}
		}
	}


	/* evolve arclength along the interface */
	bool step()
	{
		// forward Euler update
		old_nodes = nodes.dup;
		
		foreach(i, ref n; nodes)
		{
			static if(open_curve) // boundary nodes are fixed - no update.
				if(i==0 || i==nodes.length-1)
					continue;
			
			// wraparound
			size_t im1 = (i==0? nodes.length-1 : i-1);
			size_t ip1 = (i==nodes.length-1? 0 : i+1);

			// update positions
			static if(straight_springs)
			{	
				auto old_ri = old_nodes[i].r;
				auto old_rim1 = old_nodes[im1].r;
				auto old_rip1 = old_nodes[ip1].r;
				auto old_taui = old_nodes[i].tau;

				Real[2] dr_ip1, dr_i;
				dr_ip1[] = old_nodes[ip1].r[] - old_nodes[i].r[];
				dr_i[] = old_nodes[i].r[] - old_nodes[im1].r[];

				auto len_ip1 = norm(dr_ip1);
				auto len_i = norm(dr_i);

				n.s += (dt/eta)*( force_tg(len_ip1)*dotProduct(dr_ip1[],old_taui[])/len_ip1 - force_tg(len_i)*dotProduct(dr_i[], old_taui[])/len_i );

				n.s = s_PBC(n.s); // remain in [0, s_max)
				n.u = u(n.s);
			}
			else static if(curved_springs)
			{
				auto len_ip1 = s_PBC(old_nodes[ip1].s - old_nodes[i].s);
				auto len_i = s_PBC(old_nodes[i].s - old_nodes[im1].s);
				n.s += (dt/eta)*( force_tg(len_ip1) - force_tg(len_i) );
				n.s = s_PBC(n.s); // remain in [0, s_max)
				n.u = u(n.s);
			}
			else
				enforce(0, "not implemented");
			n.r[] = r(n.u);

			// update (true) tangents from position along interface
			n.tau[] = tau(n.u);

			// update length (1/density) and stress with previous node
			static if(straight_springs)
				n.length = euclideanDistance(n.r[], nodes[im1].r[]);
			else static if(curved_springs)
				n.length = s_PBC(n.s - nodes[im1].s);
			else
				enforce(0, "not implemented");
			n.stress = force_tg(n.length)/(k*a);
		}
		
		// update time
		t += dt;
		iter += 1;

		// output results
		if (iter % output_every_nth_step == 0)
		{			
			frame_iter++;
			write_state();
			stderr.writef("solving model... %d%%\r", lround(100*t/t_end));
		}
		
		return t < t_end - 0.5*dt;
	}
}

void main(string[] args)
{
	auto m = new Model;

	stderr.writeln("writing initial state...");
	m.write_state();

	stderr.write("solving model... 0%");
	while(m.step()){}
	stderr.writeln("solving model... 100%");
	
	stderr.writeln("done.");	
}
