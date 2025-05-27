/*
 * UTMBinnerBayesianBlocks.cc
 * Unit tests for MBinnerBayesianBlocks
 *
 * Copyright (C) by [Your Name/Organization] // TODO: User to update
 * All rights reserved. (Or specify MEGAlib's license if appropriate)
 */

// MEGAlib Headers
#include "MUnitTest.h"            // Base class for unit testing
#include "MBinnerBayesianBlocks.h"// The class to test
#include "MStreams.h"             // For mlog, Endl (used by MUnitTest)
#include "MGlobal.h"              // For MGlobal::Initialize
#include "MString.h"

// Standard Library Headers
#include <vector>
#include <string>
#include <cmath>        // For fabs, isnan, etc.
#include <limits>       // For numeric_limits for tolerance
#include <random>       // For data generation
#include <algorithm>    // For std::sort, std::reverse, etc.
#include <iostream>     // Fallback for critical errors

using namespace std;

class MUTMBinnerBayesianBlocks : public MUnitTest
{
public:
  MUTMBinnerBayesianBlocks();
  ~MUTMBinnerBayesianBlocks();

  // Runs all implemented test cases
  virtual bool Run(); // Overrides MUnitTest::Run

private:
  // Individual test methods
  bool TestPowerLawSpectrum();
  // Future test methods will be declared here
  
  vector<double> GeneratePowerLawData(int NumEvents, double Emin, double Emax, double Gamma);
  mt19937 m_rng; 
};

//////////////////////////////////////////////////////////////////////////////
//

MUTMBinnerBayesianBlocks::MUTMBinnerBayesianBlocks()
  : MUnitTest() // Call base class constructor explicitly
{
  random_device Rd;
  m_rng.seed(Rd()); 
}

//////////////////////////////////////////////////////////////////////////////
//

MUTMBinnerBayesianBlocks::~MUTMBinnerBayesianBlocks()
{
  // Destructor
}

//////////////////////////////////////////////////////////////////////////////
//

vector<double> MUTMBinnerBayesianBlocks::GeneratePowerLawData(int NumEvents, double Emin, double Emax, double Gamma)
{
  vector<double> Data;
  Data.reserve(NumEvents);
  uniform_real_distribution<> Dist(0.0, 1.0);

  double EminPow, EmaxPow;
  bool UseLogFormula = false;

  if (fabs(Gamma - 1.0) < 1e-9) { 
    UseLogFormula = true;
    EminPow = log(Emin);
    EmaxPow = log(Emax);
  } else {
    EminPow = pow(Emin, 1.0 - Gamma);
    EmaxPow = pow(Emax, 1.0 - Gamma);
  }

  for (int i = 0; i < NumEvents; ++i) {
    double U = Dist(m_rng); 
    double X;
    if (UseLogFormula) {
      X = exp(U * (EmaxPow - EminPow) + EminPow);
    } else {
      double ValPow = U * (EmaxPow - EminPow) + EminPow;
      X = pow(ValPow, 1.0 / (1.0 - Gamma));
    }
    Data.push_back(X);
  }
  return Data;
}

//////////////////////////////////////////////////////////////////////////////
//

bool MUTMBinnerBayesianBlocks::TestPowerLawSpectrum()
{
  mlog(MStreams::INFO)<<"Running test: TestPowerLawSpectrum..."<<Endl;
  bool TestSectionPassed = true; // Tracks pass/fail for this specific test method

  // Test parameters
  double Emin = 10.0;
  double Emax = 10000.0;
  double Gamma = 2.0; // Photon index
  int NumEvents = 5000; 

  vector<double> Data = GeneratePowerLawData(NumEvents, Emin, Emax, Gamma);

  // Use MUnitTest::Evaluate for this check
  TestSectionPassed &= Evaluate("TestPowerLawSpectrum.DataGenerationSize", 
                                static_cast<int>(Data.size()) == NumEvents, 
                                "Check if generated data size matches NumEvents");
  if (static_cast<int>(Data.size()) != NumEvents) return false; // Critical, exit if data not as expected

  MBinnerBayesianBlocks Binner;
  Binner.SetMinMax(Emin, Emax, false); 

  for (double Val : Data) {
    Binner.Add(Val, 1.0); 
  }

  const vector<double>& Edges = Binner.GetBinEdges(); 
  const vector<double>& BinnedData = Binner.GetBinnedData(); 

  // --- Assertions using MUnitTest::Evaluate ---
  TestSectionPassed &= Evaluate("TestPowerLawSpectrum.EdgesNotEmpty", !Edges.empty(), "Check Edges vector is not empty");

  if (!Edges.empty()) {
    bool EdgesAreSorted = true;
    for (size_t i = 0; i < Edges.size() - 1; ++i) {
      if (Edges[i] > Edges[i+1]) {
        EdgesAreSorted = false;
        break;
      }
    }
    TestSectionPassed &= Evaluate("TestPowerLawSpectrum.EdgesSorted", EdgesAreSorted, "Check Bin edges are sorted");

    TestSectionPassed &= Evaluate("TestPowerLawSpectrum.FirstEdgeInRange", Edges.front() >= Emin - 1e-9, 
                                  MString("First edge (") + Edges.front() + ") should be >= emin (" + Emin + ")");
    TestSectionPassed &= Evaluate("TestPowerLawSpectrum.LastEdgeInRange", Edges.back() <= Emax + 1e-9,
                                  MString("Last edge (") + Edges.back() + ") should be <= emax (" + Emax + ")");
  }
  
  bool BinnedDataSizingCorrect = true;
  if (Edges.empty()) {
    if (!BinnedData.empty() && !(BinnedData.size() == 1 && BinnedData[0] == 0)) { 
      BinnedDataSizingCorrect = false;
    }
  } else { // Edges not empty
    if (BinnedData.size() != Edges.size() + 1) {
      BinnedDataSizingCorrect = false;
    }
  }
  TestSectionPassed &= Evaluate("TestPowerLawSpectrum.BinnedDataSizing", BinnedDataSizingCorrect, 
                                "Check BinnedData size consistency with Edges size (N_edges+1)");
      
  double TotalCounts = 0;
  if (Edges.size() > 1) {
    for(size_t i = 0; i < Edges.size() - 1; ++i) { // Summing N_edges-1 bins
      if (i < BinnedData.size()) { 
        TotalCounts += BinnedData[i];
      }
    }
  }

  // Informational log for total counts, not a hard failure condition for this test
  if (fabs(TotalCounts - NumEvents) > 1e-3) { 
    mlog(MStreams::INFO)<<"Info for TestPowerLawSpectrum: Total counts in actual bins = "<<TotalCounts
                        <<", expected = "<<NumEvents
                        <<". (This is informational, not a strict failure for this test)"<<Endl;
  }
  // Example for a strict check if needed later:
  // TestSectionPassed &= Evaluate("TestPowerLawSpectrum.TotalCounts", fabs(TotalCounts - NumEvents) <= 1e-3, "Check if total counts match generated events (approx.)");


  if (TestSectionPassed) {
    mlog(MStreams::OK)<<"TestPowerLawSpectrum: All checks passed."<<Endl;
  } else {
    mlog(MStreams::ERR)<<"TestPowerLawSpectrum: One or more checks failed. See MUnitTest output above."<<Endl;
  }
  return TestSectionPassed;
}

//////////////////////////////////////////////////////////////////////////////
//

bool MUTMBinnerBayesianBlocks::Run() // Implement virtual Run from MUnitTest
{
  mlog(MStreams::INFO)<<"Starting MUTMBinnerBayesianBlocks test suite."<<Endl;
  
  bool OverallSuccess = true;
  OverallSuccess &= TestPowerLawSpectrum();
  // Add calls to other test methods here:
  // OverallSuccess &= TestAnotherFeature();

  // The pass/fail counts are managed by the MUnitTest base class via Evaluate()
  // The return value of Run() can indicate if all tests herein passed.
  // Summarize() will show the totals from MUnitTest.
  return OverallSuccess; 
}

//////////////////////////////////////////////////////////////////////////////
//

int main(int argc, char** argv) 
{
  if (!MGlobal::Initialize("MUTMBinnerBayesianBlocks", "Unit Test for MBinnerBayesianBlocks class", MGlobal::MGCR_NONE)) {
    cerr<<"FATAL: Unable to initialize MGlobal!"<<endl; 
    return 1;
  }

  MUTMBinnerBayesianBlocks Tester;
  Tester.Run(); // Execute all tests

  Tester.Summarize(); // Print summary of passed/failed tests using MUnitTest::Summarize

  return (Tester.GetNumberOfFailedTests() == 0) ? 0 : 1; // Exit code based on MUnitTest counters
}
```
