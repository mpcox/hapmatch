/*
program:	madagascar
date:		2 march 2011
author:		murray cox
function:	calculates lineages with one frequency in the ingroup, and another in the outgroup
			ingroup might be Madagascar (bottlenecked) versus outgroup Indonesia (source population)

date:		29 december 2018
update:		updated to work with libsequence2
 
compilation:
include the libsequence library
g++ -o hapmatch hapmatch.cc -lsequence -std=c++11 -Wall
*/

#include <iostream>
#include <string>
#include <vector>
#include <Sequence/SimParams.hpp>
#include <Sequence/SimDataIO.hpp> //updated for libsequence2

using namespace std;
using namespace Sequence;

void printUsage();

int main(int argc, char *argv[]){

	unsigned nruns, totsam;

	unsigned int ingroup, outgroup, thresh;
	double infreq, outfreq;

	/*
	// temporary placeholders, read in from file instead
	ingroup = 4;
	outgroup = 6;
	infreq = 0.90;
	outfreq = 0.05;
	*/

	if( argc != 6 ){
		printUsage();
		exit(1);
	}

	ingroup  = atoi(argv[1]);
	infreq   = atof(argv[2]);
	outgroup = atoi(argv[3]);
	outfreq  = atof(argv[4]);
	thresh   = atoi(argv[5]);

	/*
	// print values to screen
	cout << ingroup << "\t" << infreq << "\t" << outgroup << "\t" << outfreq << "\n";
	*/

	//define constructor to access parameters from ms program
	SimParams p;
	cin >> p;
	//cout << p << "\n";  // do not print ms command line

	nruns = p.runs();     // number of runs
	totsam = p.totsam();  // the total sample size

	// check samples sizes match user input
	if( totsam != ingroup + outgroup ){
		cout << "error: sum of ingroup and outgroup does not match simulated sample size\n";
		exit(1);
	}
    
    //SimData d( p.totsam() );  //constructor, needs to know the sample size simulated
	SimData d;  //updated for libsequence2
    
	unsigned MAXSITES = 10000;
	unsigned *indexes = static_cast< unsigned *>( malloc( MAXSITES * sizeof(unsigned) ) );

	std::ios_base::sync_with_stdio(true);

	// print out header line
	cout << "S" << "\t" << "Sin" << "\t" << "Sout" << "\t" << "hapmatch" << "\n";

	int rv;  //an integer that is the return value from fscanf, so you can check for EOF
	while( ( rv = d.fromfile(stdin) ) != EOF ){

		if( d.numsites() > MAXSITES ){  //d.numsites() --> how many positions are stored in PolyTable::positions

			MAXSITES = d.numsites() + 1;
			indexes = static_cast< unsigned *>( realloc( indexes, MAXSITES * sizeof(unsigned) )  );
		}

		/*
		// print original matrix
		cout << "original matrix:\n";
		for( unsigned int i = 0; i < totsam; ++i ){
			for( unsigned int j = 0; j < d.numsites(); ++j ){
    			cout << d[i][j];
			}
    		cout << "\n";
		}
		*/

		// reduce dataset to strings
		vector<string> all_seqs(totsam);
		for( unsigned int i = 0; i < totsam; ++i ){
			all_seqs[i] = d[i];
		}

		/*
		cout << "all_seqs:\n";
		for( unsigned int i = 0; i < totsam; ++i ){
			cout << all_seqs[i] << "\n";
		}
		*/

		// form ingroup sequences
		vector<string> ingroup_seqs(ingroup);
		for( unsigned int i = 0; i < ingroup; ++i ){
			ingroup_seqs[i] = all_seqs[i];
		}

		/*
		cout << "ingroup_seqs:\n";
		for( unsigned int i = 0; i < ingroup; ++i ){
			cout << ingroup_seqs[i] << "\n";
		}
		*/

		// form outgroup sequences
		vector<string> outgroup_seqs(outgroup);
		unsigned int seq_counter = 0;
		for( unsigned int i = ingroup; i < ingroup + outgroup; ++i ){
			outgroup_seqs[seq_counter] = all_seqs[i];
			++seq_counter;
		}

		/*
		cout << "outgroup_seqs:\n";
		for( unsigned int i = 0; i < outgroup; ++i ){
			cout << outgroup_seqs[i] << "\n";
		}
		*/

		// calculate ingroup segsites
		unsigned int ingroup_segsites = 0;
		unsigned int outgroup_segsites = 0;
		for( unsigned int a = 0; a < d.numsites(); ++a ){

			for( unsigned int x = 0; x < ingroup; ++x ){
				// use integer translation of character
				if( (int)d[x][a] == 49 ){
					++ingroup_segsites;
					break;
				}
			}

			for( unsigned int y = ingroup; y < ingroup + outgroup; ++y ){
				if( (int)d[y][a] == 49 ){
					++outgroup_segsites;
					break;
				}
			}
		}

		// generate unique ingroup strings
		vector<string> unique_ingroup_seqs(ingroup_seqs);
		std::sort(unique_ingroup_seqs.begin(), unique_ingroup_seqs.end());
		unique_ingroup_seqs.erase( unique( unique_ingroup_seqs.begin(), unique_ingroup_seqs.end() ), unique_ingroup_seqs.end() );

		vector <int>::size_type n_ingroup_haps;
		n_ingroup_haps = unique_ingroup_seqs.size();

		/*
		cout << "unique_ingroup_seqs:\n" << n_ingroup_haps << "\n";

		cout << "unique_ingroup_seqs:\n";
		for( unsigned int i = 0; i < n_ingroup_haps; ++i ){
			cout << unique_ingroup_seqs[i] << "\n";
		}
		*/

		// generate counting vector
		vector<int> c_ingroup_haps(n_ingroup_haps, 0);
		vector<int> c_outgroup_haps(n_ingroup_haps, 0);

		// count haplotypes in ingroup
		for( unsigned int i = 0; i < n_ingroup_haps; ++i ){
			for( unsigned int j = 0; j < ingroup; ++j ){

				if( unique_ingroup_seqs[i] == ingroup_seqs[j] ){
					c_ingroup_haps[i]++;
				}
			}
		}

		/*
		cout << "unique_ingroup_counts:\n";
		for( unsigned int i = 0; i < n_ingroup_haps; ++i ){
			cout << c_ingroup_haps[i] << "\n";
		}
		*/
		
		
		// generate haplotype frequencies in outgroup
		for( unsigned int i = 0; i < n_ingroup_haps; ++i ){
			if( (double)c_ingroup_haps[i]/ingroup >= infreq )
			{
				for( unsigned int j = 0; j < outgroup; ++j ){
					unsigned int x = 0;
					if( (double)c_outgroup_haps[i]/outgroup <= outfreq )
					{
						// If looking for exact matches optimise runtime by not searching char by char
						if ( thresh == 0)
						{
							if( unique_ingroup_seqs[i] == outgroup_seqs[j] )
							{
								c_outgroup_haps[i]++;
							}
						}
						
						else
						{
							for ( unsigned int q = 0; q < unique_ingroup_seqs[i].size(); ++q )
							{
								// Stop searching if threshold has been passed
								if ( x > thresh )
								{
									q = unique_ingroup_seqs[i].size();
								}
								// Search only for mismatches where the ingroup is '1' and the outgroup '0'
								if( unique_ingroup_seqs[i][q] == '1' ){
									if( outgroup_seqs[j][q] == '0'){
									++x;
									}
								}
							}
							// If the haplotypes have the desired number of mismatches count a match.
							if(	x == thresh)
							{
								c_outgroup_haps[i]++;
							}
						}
					}
					// if the maximum outgroup frequency has been reached stop checking for this ingroup sequence. 
					else {j = outgroup;}
					//cout << c_outgroup_haps[i] << "\n";
				}
			//cout << "count" << c_outgroup_haps[i] << " \n";
			}
		}

		/*
		cout << "unique_outgroup_counts:\n";
		for( unsigned int i = 0; i < n_ingroup_haps; ++i ){
			cout << c_outgroup_haps[i] << "\n";
		}

		cout << "unique_ingroup_freqs:\n";
		for( unsigned int i = 0; i < n_ingroup_haps; ++i ){
			cout << (double)c_ingroup_haps[i]/ingroup << "\n";
		}

		cout << "unique_outgroup_freqs:\n";
		for( unsigned int i = 0; i < n_ingroup_haps; ++i ){
			cout << (double)c_outgroup_haps[i]/outgroup << "\n";
		}
		*/

		// determine frequency conditions
		unsigned int matching_lineages = 0;
		for( unsigned int i = 0; i < n_ingroup_haps; ++i ){
			if( (double)c_ingroup_haps[i]/ingroup >= infreq ){
				if( (double)c_outgroup_haps[i]/outgroup <= outfreq ){
					// Only count a lineage if it is found in the outgroup as well as the ingroup.
					if (c_outgroup_haps[i] > 0)
					{
						++matching_lineages;
					}
				}
			}
		}

		// print out values
		cout << d.numsites() << "\t" << ingroup_segsites << "\t" << outgroup_segsites << "\t" << matching_lineages << "\n";

	}

	free(indexes);
}

void printUsage(){
    
    //updated: added threshold information
    cout << "usage: <ms_output> | hapmatch ingroup_size ingroup_frequency_cutoff outgroup_size outgroup_frequency_cutoff threshold"
         << "\nexample: ms 10 3 -t 4 -I 2 4 6 1 | hapmatch 4 0.444 6 0.948 0\n";
}

