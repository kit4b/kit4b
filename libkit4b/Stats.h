#pragma once
#include "commdefs.h"

/* ------------------------------------------------------------------------- 
 * This is an ANSI C library that can be used to evaluate the probability 
 * density functions (pdf's), cumulative distribution functions (cdf's), and 
 * inverse distribution functions (idf's) for a variety of discrete and 
 * continuous random variables.
 *
 * The following notational conventions are used
 *                 x : possible value of the random variable
 *                 u : real variable (probability) between 0.0 and 1.0 
 *  a, b, n, p, m, s : distribution-specific parameters
 *
 * There are pdf's, cdf's and idf's for 6 discrete random variables
 *
 *      Random Variable    Range (x)  Mean         Variance
 *
 *      Bernoulli(p)       0..1       p            p*(1-p)
 *      Binomial(n, p)     0..n       n*p          n*p*(1-p)
 *      Equilikely(a, b)   a..b       (a+b)/2      ((b-a+1)*(b-a+1)-1)/12 
 *      Geometric(p)       0...       p/(1-p)      p/((1-p)*(1-p))
 *      Pascal(n, p)       0...       n*p/(1-p)    n*p/((1-p)*(1-p))
 *      Poisson(m)         0...       m            m
 *
 * and for 7 continuous random variables
 *
 *      Uniform(a, b)      a < x < b  (a+b)/2      (b-a)*(b-a)/12
 *      Exponential(m)     x > 0      m            m*m
 *      Erlang(n, b)       x > 0      n*b          n*b*b
 *      Normal(m, s)       all x      m            s*s
 *      Lognormal(a, b)    x > 0         see below
 *      Chisquare(n)       x > 0      n            2*n
 *      Student(n)         all x      0  (n > 1)   n/(n-2)   (n > 2)
 *
 * For the Lognormal(a, b), the mean and variance are
 *
 *                        mean = Exp(a + 0.5*b*b)
 *                    variance = (Exp(b*b) - 1)*Exp(2*a + b*b)
 *
 * Name            : rvms.c (Random Variable ModelS)
 * Author          : Steve Park & Dave Geyer
 * Language        : ANSI C
 * Latest Revision : 11-22-97
 * ------------------------------------------------------------------------- 
 */

// Wald–Wolfowitz runs test table
// Required - number of runs 'r', number of values 'n1' and number of values 'n2'
// critical value is the minimum number of runs at a Pvalue >= 0.025
const uint8_t WaldWolfowitzRunsCritValues[19][19] = {
{1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2},
{1,1,1,1,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3},
{1,1,1,2,2,2,3,3,3,3,3,3,3,3,4,4,4,4,4},
{1,1,2,2,3,3,3,3,3,4,4,4,4,4,4,4,5,5,5},
{1,2,2,3,3,3,3,4,4,4,4,5,5,5,5,5,5,6,6},
{1,2,2,3,3,3,4,4,5,5,5,5,5,6,6,6,6,6,6},
{1,2,3,3,3,4,4,5,5,5,6,6,6,6,6,7,7,7,7},
{1,2,3,3,4,4,5,5,5,6,6,6,7,7,7,7,8,8,8},
{1,2,3,3,4,5,5,5,6,6,7,7,7,7,8,8,8,8,9},
{1,2,3,4,4,5,5,6,6,7,7,7,8,8,8,9,9,9,9},
{2,2,3,4,4,5,6,6,7,7,7,8,8,8,9,9,9,10,10},
{2,2,3,4,5,5,6,6,7,7,8,8,9,9,9,10,10,10,10},
{2,2,3,4,5,5,6,7,7,8,8,9,9,9,10,10,10,11,11},
{2,3,3,4,5,6,6,7,7,8,8,9,9,10,10,11,11,11,12},
{2,3,4,4,5,6,6,7,8,8,9,9,10,10,11,11,11,12,12},
{2,3,4,4,5,6,7,7,8,9,9,10,10,11,11,11,12,12,13},
{2,3,4,5,5,6,7,8,8,9,9,10,10,11,11,12,12,13,13},
{2,3,4,5,6,6,7,8,8,9,10,10,11,11,12,12,13,13,13},
{2,3,4,5,6,6,7,8,9,9,10,10,11,12,12,13,13,13,14}
};



const int cMaxChiSqrRows = 10;	  // max number of rows accepted by CalcChiSqr
const int cMaxChiSqrCols = 250;	  // max number of columns accepted by CalcChiSqr
const int cAllocLogFacts = 500000;  // additional chunk alloc size when allocating for m_pLogFact[]
const int cMaxTotSampleCnt = 150000000;	// scale down rows such that the total sample count over all rows is less than this

const double SQRT2PI = 2.506628274631;		/* sqrt(2 * pi) */
const double TINY =	1.0e-10;

class CStats
{


	double m_lncof[6];
	int m_NumLogFacts;		// largest precalculated log(factorial) in m_pLogFact
	int m_AllocLogFacts;	// how many entries have been allocated in m_pLogFact
	double *m_pLogFact;		// holds array of all factorials from 1..m_AllocLogFacts as log(fact)
	double lngam(double z);
	double gamminc(double aa, double xx);



public:
	CStats(void);
	~CStats(void);
	void Init(void);
	double					// returned P1 
		FishersExactTest(int R1C1,		// sample1 true
						 int R1C2,		// sample1 false
						 int R2C1,		// sample2 true
						 int R2C2);		// sample2 false;

	double ChiSqr2PVal(int df, double ChiSqr); // returns P-value

	double							// returned Chi-Square, -1.0 if any expected count is less than 5
		CalcChiSqr(int Rows,		// number of rows in pCells table
				   int Cols,		// number of columns in pCells table
				   int *pCells);	// expected to contain [Rows][Cols] counts - e.g. values for column 1, followed by values for column 2...
	
//	double Rand(void);				// returns random value 0.0 <= Rand < 1.0

	double Calc_nCk(uint32_t n, uint32_t k); // Calculates nCk = n! / (n-k)!k!

	double ProbKeqlk(uint32_t n, uint32_t k, double p); //  Calculates Pr(K = k) = nCk * p^k * q^(n-k)

	// The cumulative distribution is the sum of all the success probabilities K = 0 up to K = k in n trials
	// Note: left tailed
	double Binomial(int n, // number of trials 
				int k,	// number of observed successes 
				double p); // expected probability of a success in a trial

	bool	// true if single sided Wald–Wolfowitz runs test shows support at >= 0.025 for sequence of called haplotypes randomly drawn from FaFb
			IsRandomHaplotypesFaFb(uint32_t n1,	// number of Fa haplotype calls
				  uint32_t n2,	// number of Fb haplotype calls
				  uint32_t r);	// number of runs in sequence of haplotypes

	double phi(double x); // one tailed quick and dirty cumulative density function (CDF) of a standard normal (Gaussian) random variable.


	double LogFactorial(long n);
	double LogChoose(long n, long m);

	double pdfBernoulli(double p, long x);
	double cdfBernoulli(double p, long x);
	long   idfBernoulli(double p, double u);

	double pdfEquilikely(long a, long b, long x);
	double cdfEquilikely(long a, long b, long x);
	long   idfEquilikely(long a, long b, double u);

	double pdfBinomial(long n, double p, long x);
	double cdfBinomial(long n, double p, long x);
	long   idfBinomial(long n, double p, double u);

	double pdfGeometric(double p, long x);
	double cdfGeometric(double p, long x);
	long   idfGeometric(double p, double u);

	double pdfPascal(long n, double p, long x);
	double cdfPascal(long n, double p, long x);
	long   idfPascal(long n, double p, double u);

	double pdfPoisson(double m, long x);
	double cdfPoisson(double m, long x);
	long   idfPoisson(double m, double u);

	double pdfUniform(double a, double b, double x);
	double cdfUniform(double a, double b, double x);
	double idfUniform(double a, double b, double u);

	double pdfExponential(double m, double x);
	double cdfExponential(double m, double x);
	double idfExponential(double m, double u);

	double pdfErlang(long n, double b, double x);
	double cdfErlang(long n, double b, double x);
	double idfErlang(long n, double b, double u);

	double pdfNormal(double m, double s, double x);
	double cdfNormal(double m, double s, double x);
	double idfNormal(double m, double s, double u);

	double pdfLognormal(double a, double b, double x);
	double cdfLognormal(double a, double b, double x);
	double idfLognormal(double a, double b, double u);

	double pdfChisquare(long n, double x);
	double cdfChisquare(long n, double x);
	double idfChisquare(long n, double u);

	double pdfStudent(long n, double x);
	double cdfStudent(long n, double x);
	double idfStudent(long n, double u);

	double pdfStandard(double x);
	double cdfStandard(double x);
	double idfStandard(double u);

	double LogGamma(double a);
	double LogBeta(double a, double b);
	double InBeta(double a, double b, double x);
	double InGamma(double a, double x);
};
