/*
 * MUnitTest.h
 *
 * Copyright (C) by Andreas Zoglauer.
 * All rights reserved.
 *
 * Please see the source-file for the copyright-notice.
 *
 */


#ifndef __MUnitTest__
#define __MUnitTest__


////////////////////////////////////////////////////////////////////////////////


// Standard libs:
#include <vector> 

// MEGAlib libs:
#include "MGlobal.h"
#include "MStreams.h" 
#include "MString.h"  

using namespace std; 


////////////////////////////////////////////////////////////////////////////////


//! The base class for unit test
class MUnitTest
{
public:
  MUnitTest();
  virtual ~MUnitTest();
  
  template <typename T1, typename T2> bool Evaluate(MString Function, T1 Input, MString Description, T2 Output, T2 Truth)
  {
    if (Output != Truth) {
      mlog(MStreams::ERR)<<"FAILED: "<<Function<<"  <-- "<<Input<<Endl;
      mlog(MStreams::ERR)<<"   Description: "<<Description<<Endl;
      mlog(MStreams::ERR)<<"   Expected:    "<<Truth<<Endl;
      mlog(MStreams::ERR)<<"   Output:      "<<Output<<Endl;
      ++m_NumberOfFailedTests;
      return false;
    }
    ++m_NumberOfPassedTests;
    return true;
  }

  // Declaration only:
  bool Evaluate(MString TestName, bool Condition, MString FailureDescription);

  // Declaration only:
  bool Evaluate(MString TestName, const vector<double>& Output, const vector<double>& Truth, MString Description, double Tolerance = 1e-6);

  virtual bool Run() = 0;
  void Summarize();

  unsigned int GetNumberOfFailedTests() const { return m_NumberOfFailedTests; }
  unsigned int GetNumberOfPassedTests() const { return m_NumberOfPassedTests; }
  
private:
   unsigned int m_NumberOfPassedTests;
   unsigned int m_NumberOfFailedTests;

#ifdef ___CLING___
 public:
  ClassDef(MUnitTest, 1) 
#endif
};
#endif
