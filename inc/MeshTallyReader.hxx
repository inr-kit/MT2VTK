#ifndef MESHTALLYREADER_HXX
#define MESHTALLYREADER_HXX


// // //
//       EXTERNAL INCLUDES
/////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include <fstream>
#include <map>
#include <list>

#include <McnpMeshGeometryType.hxx>
#include <McnpTally.hxx>

using namespace std;


// // //
//       DEFINITIONS
//////////////////////////////////////////////


enum WriteType{
    Write_DATA,
    Write_ERROR
};

enum IntegrateConst{
    CONST_X,
    CONST_Y,
    CONST_Z
};


/*!
  MeshTallyReader
  ===============
  Purpose: The MeshTallyReader reads a MeshTally file produced by MCNP and processes the containing data
           The processed data can be translated to a 3D CAD model with colormaps which correspond to the
           tally values (heatmap)
*/
class MeshTallyReader {
public:
  // ctor and dtor
    MeshTallyReader();
    // MeshTallyReader(const MeshTallyReader& orig);
    MeshTallyReader(string& inFileName);
    virtual ~MeshTallyReader();

 // File operations
//! sets the name of the mesh tally file produced by MCNP
   void SetInFileName(string& fileName);

//! returns the file name
   string GetInFileName();

//! Read and process Tally File
   bool ReadTallyFile(string& fileName);

//! Read and process Tally File - myInFileName must be set previously
   bool ReadTallyFile();

//! Set the file name of the output STEP file
   void SetOutFileName(string& fileName);

//! returns the file name of the ouput STEP file
   string GetOutFileName();

//! Write data to VTK file - myOutFileName must be defined
   bool WriteVTK();

//! Write data to VTK file
   bool WriteVTK(string& outName);

//! Write selection of tallies + errors
   bool WriteSelectedTallies(string& selectionList);

//! Write All Tallies
   bool WriteAll();

//! Write relative errors to VTK file
   bool WriteErrorVTK();

//! Write relative errors to VTK file
   bool WriteErrorVTK(string& outName);

//!  Multiplies each value of the given mesh with a constant factor
   void Multiply(int meshNumber, float factor);

//! call Multiply(int,int) for all bound meshes
   void Multiply(float factor);

//! Returns number of all meshes in all tallies
   int NumberOfMeshes() const;

//! Returns Number Of Tallies found in file
   int NumberOfTallies() const;

//! performs a mathematical operation between two tallies
   // e.g. 14 + 4 would add the results of mesh 14 with results of mesh 4
   //      14 / 4 would divide the results of 14 by the results of 4
   void PerformOperation(string strTally1, string operation, string strTally2=string("0"));

//! Print Tally Information to screen
   void PrintTallyInformation() const;

//! Print Available Meshes (Tallies and Energy Bins)
   void PrintMeshesAvailable() const;

//! Print Available Tallies
   void PrintTalliesAvailable();

//! Integrates all values of constant x or y or z
   void MakeIntegralList();

//! Adding a Tally to the tally map
   bool AddTally(McnpTally& tally);

//! Getting a Tally by TallyNumber
   McnpTally GetTally(int tallyNumber);


//! Set scaling factor for geometry : useful if original units were mm but calculation results (meshtal) are in cm
   void SetGeometryScalingFactor(const float& scalingFactor);

//! Get scaling factor for geometry
   float GetGeometryScalingFactor() const;

// STATIC FUNCTIONS
////////////////////

   static void PerformOperationOnTallies(McnpTally* tal1, string operation, McnpTally* tal2);




private:

// FUNCTIONS
/////////////////
   void Init();
   string GetLine(ifstream& fileStream);
   int GetTallyNumber(ifstream& fileStream);
   McnpParticleType GetParticleType(ifstream& fileStream);
   void GetBinBoundaries(ifstream& fileStream, McnpTally& curTalDat);
   bool StringToFloat(string str, float& val);
   void GetEnergyBinBoundaries(ifstream& fileStream, McnpTally& curTalDat);
   McnpMeshFormat DetermineTallyFormat(ifstream& fileStream, McnpTally& curTal);
   void GetColoumnData(ifstream& fileStream, McnpTally& curTal);
   void GetMatrixData(ifstream& fileStream, McnpTally& curTal);
   void PrintProgress(const int& n, const int& N) const;
   bool Write(int tallyNumber, WriteType wType);
   bool CheckTallies();

// FIELDS
/////////////////
   map<int, McnpTally> myTallyMap;

   string myInFileName;
   string myOutFileName;
   int myLineCount;
   float myGeometryScalingFactor;
};

#endif // MESHTALLYREADER_HXX
