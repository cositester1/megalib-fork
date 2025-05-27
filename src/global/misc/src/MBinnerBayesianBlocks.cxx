/*
 * MBinnerBayesianBlocks.cxx
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
#include "MBinnerBayesianBlocks.h"

// Standard libs:
#include <cmath>
#include <algorithm>
#include <sstream> // For using ostringstream

// ROOT libs:

// MEGAlib libs:
#include "MGlobal.h"


////////////////////////////////////////////////////////////////////////////////


#ifdef ___CLING___
ClassImp(MBinnerBayesianBlocks)
#endif


////////////////////////////////////////////////////////////////////////////////


//! Default constructor
MBinnerBayesianBlocks::MBinnerBayesianBlocks() : m_MinimumBinWidth(0.000001), m_MinimumCountsPerBin(0), m_Prior(4), m_UseBinning(true)
{
}


////////////////////////////////////////////////////////////////////////////////


//! Default destructor
MBinnerBayesianBlocks::~MBinnerBayesianBlocks()
{
}


////////////////////////////////////////////////////////////////////////////////


void Print(vector<double>& Array) {
  std::ostringstream oss;
  for (unsigned int i = 0; i < Array.size(); ++i) {
    oss << Array[i] << " ";
  }
  MGlobal::msg(MGlobal::DEB) << "MBinnerBayesianBlocks::Print(double): " << oss.str() << Gendl;
}


////////////////////////////////////////////////////////////////////////////////


void Print(vector<int>& Array) {
  std::ostringstream oss;
  for (unsigned int i = 0; i < Array.size(); ++i) {
    oss << Array[i] << " ";
  }
  MGlobal::msg(MGlobal::DEB) << "MBinnerBayesianBlocks::Print(int): " << oss.str() << Gendl;
}


////////////////////////////////////////////////////////////////////////////////


//! The actual histogramming process - default just makes one bin
void MBinnerBayesianBlocks::Histogram()
{
  if (m_IsModified == false) return;
  
  m_BinEdges.clear();
  m_BinnedData.clear();
  
  // Step 1: Sort the Array increasing
  m_Values.sort(SortBinnedData);

  double Min = m_Minimum;
  double Max = m_Maximum;
  if (m_Adapt == true) {
    double Front = m_Values.front().m_AxisValue;
    double Back = m_Values.back().m_AxisValue;
    if (Front < Back) {
      if (Front > Min && Front < Max) Min = Front;
      if (Back < Max && Back > Min) Max = Back;
    }
  }
    
  unsigned int Size = 0;

  // Step 2: Create cell edges
  vector<double> Edges;
  Edges.push_back(Min);
  if (m_UseBinning == true) {
    while (Edges.back() < Max) {
      Edges.push_back(Edges.back() + m_MinimumBinWidth);
    }
    if (Edges.back() < Max) Edges.push_back(Max);
    Size = Edges.size() - 1;
  } else {
    Size = m_Values.size();
    MBinnedData Last = m_Values.front();
    for (list<MBinnedData>::iterator I = ++(m_Values.begin()); I != m_Values.end(); ++I) {
      Edges.push_back(0.5*(Last.m_AxisValue + (*I).m_AxisValue));
      Last = (*I);
    }
    Edges.push_back(Max);
  }
  //cout<<"Edges:"<<endl;
  //Print(Edges);

  // Step 3: Create Block length:
  vector<float> BlockLength;
  for (unsigned int i = 0; i < Edges.size(); ++i) {
    BlockLength.push_back(Edges.back() - Edges[i]);
  }
  //cout<<"Block length:"<<endl;
  //Print(BlockLength);
  
  // Step 4: Prepare for iterations
  vector<float> CountsPerBin(Size, 0);
  for (list<MBinnedData>::iterator I = m_Values.begin(); I != m_Values.end(); ++I) {
    double Value = (*I).m_AxisValue;
    for (unsigned int e = 0; e < Edges.size() - 1; ++e) { // Speed improvement possible: Since m_Values and Edges are sorted, a linear scan (merge-like) approach (O(N_values + N_edges)) can be used instead of this O(N_values * N_edges) nested loop.
      if (Edges[e] <= Value && Edges[e+1] > Value) {
        CountsPerBin[e] += (*I).m_DataValue;
        break;     
      }
    }
  }
  //cout<<"Counts:"<<endl;
  //Print(CountsPerBin);
  
  vector<float> Best(Size, 0.0);
  vector<unsigned int> Last(Size, 0);

  // Step 5: Iterate
  for (unsigned int s = 0; s < Size; ++s) {
    //cout<<s<<" / "<<Size<<endl;
  
    // Calculate the width of the blocks
    vector<float> Width; // log(float) is the fastest of the log calculations
    for (unsigned int i = 0; i <= s; ++i) {
      Width.push_back(BlockLength[i] - BlockLength[s+1]);
    }
    //cout<<"Width: "<<endl;
    //Print(Width);
    
    // Calculate the block count
    vector<float> BlockCounts(s + 1, 0);
    float current_sum = 0;
    for (int k = s; k >= 0; --k) { // Iterate downwards from s to 0
      current_sum += CountsPerBin[k];
      BlockCounts[k] = current_sum;
    }
    //cout<<"BlockCounts: "<<endl;
    //Print(BlockCounts);

    //
    vector<float> Fits;
    for (unsigned int i = 0; i <= s; ++i) {
      float Fit;
      if (Width[i] <= 0) {
        Fit = -1.0e38f; // A very large negative number
      } else {
        if (BlockCounts[i] == 0) {
          // N log N term is 0. The formula is N (log N - log W) = N log N - N log W.
          // So if N=0, this becomes 0 - 0 * log W = 0.
          Fit = 0.0f;
        } else {
          // Both BlockCounts[i] and Width[i] are positive
          Fit = BlockCounts[i] * (logf(BlockCounts[i]) - logf(Width[i]));
        }
      }
      Fit -= m_Prior;
      Fits.push_back(Fit);
    }
    //cout<<"Fits (2): "<<endl;
    //Print(Fits);
    for (unsigned int i = 1; i <= s; ++i) {
      Fits[i] += Best[i-1];
    }
    //cout<<"Fits (3): "<<endl;
    //Print(Fits);
  
    unsigned int Maximum = 0;
    for (unsigned int i = 0; i < Fits.size(); ++i) {
      if (Fits[i] > Fits[Maximum]) Maximum = i;
    }
    Last[s] = Maximum;
    Best[s] = Fits[Maximum];
  }
  
  // Scargle's implementation breaks when Size == 1
  // Step 6: Find the change points
  std::vector<unsigned int> GeneratedChangePoints;
  unsigned int CurrentCPIndex = Size; // 'Size' is the number of cells, cell indices 0 to Size-1. Edges indices 0 to Size.

  while (CurrentCPIndex > 0) {
      GeneratedChangePoints.push_back(CurrentCPIndex);

      // 'Last' has 'Size' elements, valid indices 0 to Size-1.
      // CurrentCPIndex ranges from 'Size' down to 1.
      // So, CurrentCPIndex - 1 ranges from 'Size - 1' down to 0, which is a valid index for 'Last'.
      unsigned int PrevCPIndex = Last[CurrentCPIndex - 1];

      if (PrevCPIndex >= CurrentCPIndex) {
          MGlobal::msg(MGlobal::ERR) << "Bayesian Blocks change point reconstruction failed. Index did not decrease: Current=" << CurrentCPIndex << ", Previous=" << PrevCPIndex << Gendl;
          GeneratedChangePoints.clear(); // Invalidate results
          break;
      }
      CurrentCPIndex = PrevCPIndex;
  }

  if (CurrentCPIndex == 0) {
      GeneratedChangePoints.push_back(0); // Add the starting point (index 0 for Edges array)
  } else {
      // This block is reached if the loop broke prematurely due to an error (e.g., PrevCPIndex >= CurrentCPIndex)
      // or if Size was initially 0 and the loop didn't run (CurrentCPIndex remains 0, so this else is skipped).
      if (!GeneratedChangePoints.empty()) { // Implies an error occurred and we broke from the loop
          MGlobal::msg(MGlobal::ERR) << "Something went wrong with the change points during Bayesian Block binning. Path to 0 not found." << Gendl;
      }
      // If GeneratedChangePoints is empty here AND Size > 0, it means the error happened on the first try.
      // If Size was 0, CurrentCPIndex starts as 0, loop is skipped, CurrentCPIndex == 0 is true, {0} is added.
  }

  std::reverse(GeneratedChangePoints.begin(), GeneratedChangePoints.end());

  m_BinEdges.clear();

  if (GeneratedChangePoints.empty() && Size > 0) {
      MGlobal::msg(MGlobal::WAR) << "No change points were generated by Bayesian Blocks for Size = " << Size << ". Resulting binning may be trivial." << Gendl;
      // As a fallback, one might consider adding Edges[0] and Edges[Size] if Edges is not empty and Size > 0
      // if (Edges.size() > Size && Size > 0) { // Edges should have Size+1 elements
      //    m_BinEdges.push_back(Edges[0]);
      //    m_BinEdges.push_back(Edges[Size]);
      // }
  }

  for (unsigned int PointIndex : GeneratedChangePoints) {
      if (PointIndex < Edges.size()) { // Edges contains original cell boundaries. PointIndex is an index for Edges.
          m_BinEdges.push_back(Edges[PointIndex]);
      } else {
          MGlobal::msg(MGlobal::ERR) << "Change point index " << PointIndex << " is out of bounds for Edges vector (size " << Edges.size() << ")." << Gendl;
          m_BinEdges.clear(); // Critical error, invalidate m_BinEdges
          break;
      }
  }
  
  // Step 7: Do some sanity checks:
  if (m_UseBinning == false) {
    // Reject bins which are smaller than X
    if (m_BinEdges.size() > 2) {
      for (unsigned int i = 1; i < m_BinEdges.size(); ++i) {
        //cout<<m_BinEdges[i-1]<<".."<<m_BinEdges[i]<<endl;
        if (m_BinEdges[i] - m_BinEdges[i-1] < m_MinimumBinWidth) {
          if (i == 1) {
            m_BinEdges.erase(m_BinEdges.begin()+i);
            i--;
          } else if (i == m_BinEdges.size() - 1) {
            m_BinEdges.erase(m_BinEdges.begin() + i - 1); // Corrected to use m_BinEdges.begin() for clarity, same as m_BinEdges.end()-2
            i--;
          } else { // An intermediate bin is too small, merge with right neighbor
            m_BinEdges.erase(m_BinEdges.begin()+i);
            i--;
          }
        }
      }
      //for (unsigned int i = 1; i < m_BinEdges.size(); ++i) {
      //  cout<<m_BinEdges[i]<<endl;
      //}
    }
  }
    
  // Step 8: Finally fill the data array
  m_BinnedData.resize(m_BinEdges.size()+1, 0);
  for (list<MBinnedData>::iterator I = m_Values.begin(); I != m_Values.end(); ++I) {
    for (unsigned int e = 0; e < m_BinEdges.size(); ++e) {
      if (m_BinEdges[e] > (*I).m_AxisValue) {
        if (e > 0) {
          m_BinnedData[e-1] += (*I).m_DataValue;
        }
        break;
      }
    }
  }    

  // Step 9: Reject bins with less than X elements
  if (m_MinimumCountsPerBin > 0) {
    for (unsigned int e = 0; e < m_BinEdges.size()-1; ++e) {
      //cout<<"Content: "<<m_BinnedData[e]<<" going from "<<m_BinEdges[e]<<" - "<<m_BinEdges[e+1]<<endl;
      if (m_BinnedData[e] < m_MinimumCountsPerBin && e < m_BinEdges.size() - 1) {
        //cout<<"Erasing..."<<endl;
        // Move higher content down and erase bins
        m_BinnedData[e] += m_BinnedData[e+1];
        m_BinnedData.erase(m_BinnedData.begin()+e+1);
        m_BinEdges.erase(m_BinEdges.begin()+e+1);
        e--;
      }
    }
  }
  
  
  m_IsModified = false;
}


// MBinnerBayesianBlocks.cxx: the end...
////////////////////////////////////////////////////////////////////////////////
