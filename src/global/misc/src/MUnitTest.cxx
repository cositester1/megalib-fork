/*
 * MUnitTest.cxx
 *
 *
 * Copyright (C) by Andreas Zoglauer.
 * All rights reserved.
 *
 *
 * This code implementation is the intellectual property of
 * Andreas Zoglauer.
 *
 * By copying, distributing or modifying the Program (or any work
 * based on the Program) you indicate your acceptance of this statement,
 * and all its terms.
 *
 */


// Include the header:
#include "MUnitTest.h"
#include <cmath> // For fabs in vector<double> Evaluate

// using namespace std; is in MUnitTest.h
// MStreams.h and MGlobal.h are included via MUnitTest.h

////////////////////////////////////////////////////////////////////////////////
#ifdef ___CLING___
ClassImp(MUnitTest)
#endif
////////////////////////////////////////////////////////////////////////////////

MUnitTest::MUnitTest() :
  m_NumberOfPassedTests(0),
  m_NumberOfFailedTests(0)
{
  // Constructor body can be empty if all done in initializer list
}

////////////////////////////////////////////////////////////////////////////////

MUnitTest::~MUnitTest()
{
  // Nothing to do here
}

//////////////////////////////////////////////////////////////////////////////
//

bool MUnitTest::Evaluate(MString TestName, bool Condition, MString FailureDescription)
{
  if (!Condition) {
    mlog(MStreams::ERR)<<"FAILED: "<<TestName<<Endl;
    mlog(MStreams::ERR)<<"   Condition not met: "<<FailureDescription<<Endl;
    ++m_NumberOfFailedTests;
    return false;
  }
  ++m_NumberOfPassedTests;
  return true;
}

//////////////////////////////////////////////////////////////////////////////
//

bool MUnitTest::Evaluate(MString TestName, const vector<double>& Output, const vector<double>& Truth, MString Description, double Tolerance)
{
  bool Match = true;
  if (Output.size() != Truth.size()) {
    Match = false;
  } else {
    for (size_t i = 0; i < Output.size(); ++i) {
      if (fabs(Output[i] - Truth[i]) > Tolerance) {
        Match = false;
        break;
      }
    }
  }

  if (!Match) {
    mlog(MStreams::ERR)<<"FAILED: "<<TestName<<" - Vector comparison failure."<<Endl;
    mlog(MStreams::ERR)<<"   Description: "<<Description<<Endl;
    // Could add detailed vector logging here if desired for very verbose failures
    ++m_NumberOfFailedTests;
    return false;
  }
  ++m_NumberOfPassedTests;
  return true;
}

//////////////////////////////////////////////////////////////////////////////
//

//! Summarize the test run
void MUnitTest::Summarize()
{
  mlog(MStreams::INFO)<<"--- Test Summary ---"<<Endl;
  mlog(MStreams::INFO)<<"Passed tests: "<<m_NumberOfPassedTests<<Endl;
  mlog(MStreams::INFO)<<"Failed tests: "<<m_NumberOfFailedTests<<Endl;
  
  if (m_NumberOfFailedTests == 0 && m_NumberOfPassedTests > 0) {
    mlog(MStreams::OK)<<"All "<<m_NumberOfPassedTests<<" tests passed successfully."<<Endl;
  } else if (m_NumberOfPassedTests == 0 && m_NumberOfFailedTests == 0) {
    mlog(MStreams::WARN)<<"No tests were run or recorded."<<Endl;
  } else {
    mlog(MStreams::ERR)<<"There were "<<m_NumberOfFailedTests<<" failed tests out of "<< (m_NumberOfPassedTests + m_NumberOfFailedTests) << " total tests."<<Endl;
  }
  mlog(MStreams::INFO)<<"--------------------"<<Endl;
}


// MUnitTest.cxx: the end...
////////////////////////////////////////////////////////////////////////////////
