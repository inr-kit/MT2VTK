#include "MeshTallyReader.hxx"

// std Classes
//////////////////////////////
#include <sstream>
#include <fstream>
#include <vector>
#include <iomanip>

// boost Classes
////////////////////////////////
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <boost/math/special_functions/trunc.hpp>
#include <boost/math/special_functions/round.hpp>
#include <boost/algorithm/string.hpp>
typedef boost::tokenizer<boost::char_separator<char> > tokenizer;

bool MeshTallyReader::CheckTallies()
{
    for(map<int, McnpTally>::iterator it = myTallyMap.begin(); it != myTallyMap.end(); it++) // Tallies
    {
        McnpTally curTal = (*it).second;
        for(int e=0; e<curTal.EnergyDimension()-1; e++)
            for(int i=0; i<curTal.XDimension(); i++)
                for(int j=0; j<curTal.YDimension(); j++)
                    for(int k=0; k<curTal.ZDimension(); k++)
                        if(curTal.GetValue(e,i,j,k) < 0 || float (fabs(curTal.GetValue(e,i,j,k)) > 1.e30))
                        {
                            cout << "\nFailure for : " << curTal.GetTallyNumber() <<
                                    "->(" << e <<") [" << i << ", " << j << ", " <<
                                    k <<"] = " << curTal.GetValue(e,i,j,k) << endl;
                            return false;
                        }
    }
    return true;
}


/////////////////////////////////////
//
//   CTORS AND DTORS
//
/////////////////////////////////////

MeshTallyReader::MeshTallyReader()
{
    Init();
}

MeshTallyReader::MeshTallyReader(string& inFileName)
{
    Init();
    myInFileName = inFileName;
    myLineCount = 0;
    if(!ReadTallyFile())
        cout << "ERROR :: Tally File: " << myInFileName << "could not be read/evaluated!!!\n";
}

/*MeshTallyReader::MeshTallyReader(const MeshTallyReader& orig)
{
    Init();
}*/


MeshTallyReader::~MeshTallyReader()
{

}

/////////////////////////////////////////////////////////////


//////////////////////////////////////////////////
//
//  FUNCTIONS
//
//////////////////////////////////////////////////

/// ADD TALLY
/// //////////////////////
bool MeshTallyReader::AddTally(McnpTally &tally)
{
    if(tally.GetTallyNumber() < 0)
        return false;

    myTallyMap[tally.GetTallyNumber()] = tally;

    return true;
}

/// GET INFILE NAME
/// //////////////////////
string MeshTallyReader::GetInFileName()
{
    return myInFileName;
}

/// GET OUTFILE NAME
/// //////////////////////
string MeshTallyReader::GetOutFileName()
{
    return myOutFileName;
}



/// MAKE INTEGRAL LIST -- too specialised for general purposes!!!
/// ///////////////////////
void MeshTallyReader::MakeIntegralList()
{
    list<float> integralXValuesLowerBin;
    list<float> integralXValuesUpperBin;
    list<float> integralXValuesTotal;

    int xCnt(0), yCnt(0), zCnt(0);

    ofstream printFile("integralYZ.dat");
    cout << "Writing x integral values to : integralYZ.dat ...\n";

    for(map<int,McnpTally>::const_iterator it = myTallyMap.begin(); it != myTallyMap.end(); it++ )
    {
        McnpTally tal = (*it).second;

        list<float>::const_iterator itX = tal.GetXBinBegin();

        // RESTRICTION BLOCK
       /* bool restricted(false);
        float xMin(830.0f), xMax(1500.0f), xTurningPoint(1178.0f);
        float yMin(-158.0f), yMax(158.0f);
        float zMin(-121.0f), zMax(245.0f);
        float yMinRest(-107.0f), yMaxRest(107.0f);
        float zMinRest(-70.0f), zMaxRest(190.0f);*/

        for( xCnt = 0; xCnt < tal.XDimension() ; xCnt++ , itX++)
        {
            float sumLowerBinYZ(0.0f), sumUpperBinYZ(0.0f), sumTotalYZ(0.0f);

            // RESTRICTION BLOCK
           /* if( (*itX) < xMin || (*itX) > xMax )
                continue;
            if( (*itX) < xTurningPoint )
                restricted = true;
            else
                restricted = false;*/
            list<float>::const_iterator itY = tal.GetYBinBegin();

            for(yCnt = 0; yCnt < tal.YDimension(); yCnt++, itY++ )
            {
                // RESTRICTION BLOCK
                /*if( restricted )
                {
                    if( (*itY) < yMinRest || (*itY) > yMaxRest )
                        continue;
                }
                else if( (*itY) < yMin || (*itY) > yMax )
                    continue;*/
                list<float>::const_iterator itZ = tal.GetZBinBegin();

                for(zCnt = 0; zCnt < tal.ZDimension(); zCnt++, itZ++ )
                {
                    //RESTRICTION BLOCK
                    /*if( restricted )
                    {
                        if( (*itZ) < zMinRest || (*itZ) > zMaxRest )
                            continue;
                    }
                    else if( (*itZ) < zMin || (*itZ) > zMax )
                        continue;*/

                    sumLowerBinYZ += tal.GetValue(0, xCnt, yCnt, zCnt);
                    sumUpperBinYZ += tal.GetValue(1, xCnt, yCnt, zCnt);
                    sumTotalYZ += tal.GetValue(2, xCnt, yCnt, zCnt);
                }
            }

            printFile << "#\n# " << myInFileName.c_str() << "\n#\n";
            printFile << "# X \t Lower Energy\t Upper Energy\t Total\n";
            printFile << "# ---------------------------------------------------------------\n";
            printFile << (*itX) << "\t" << sumLowerBinYZ << "\t" << sumUpperBinYZ << "\t" << sumTotalYZ << endl;

            integralXValuesLowerBin.push_back(sumLowerBinYZ);
            integralXValuesUpperBin.push_back(sumUpperBinYZ);
            integralXValuesTotal.push_back(sumTotalYZ);
        }
        break;
    }
}


/// MULTIPLY
/// //////////////////////
void MeshTallyReader::Multiply(int meshNumber, float factor)
{
    McnpTally curTal = myTallyMap[meshNumber];
    curTal.Multiply(factor);
}

/// MULTIPLY
/// //////////////////////
void MeshTallyReader::Multiply(float factor)
{
    if(NumberOfTallies() <= 0)
    {
        cout << "No tallies bound in tally map. Make sure to call ReadTallyFile() first\n";
        return;
    }

    for(map<int, McnpTally>::iterator it = myTallyMap.begin(); it != myTallyMap.end(); it++)    
        Multiply((it->second).GetTallyNumber(), factor);    
}

///  NUMBER OF MESHES
/// ///////////////////////
int MeshTallyReader::NumberOfMeshes() const
{
    int numberOfMeshes(0);

    for(map<int,McnpTally>::const_iterator it=myTallyMap.begin(); it!=myTallyMap.end(); it++)
    {

        numberOfMeshes += ( ((*it).second).EnergyDimension()-1 );
    }

    cout << "NUMBER OF MESHES : " << numberOfMeshes << endl;

    return numberOfMeshes;
}


///  NUMBER OF Tallies
/// //////////////////////
int MeshTallyReader::NumberOfTallies() const
{
    return myTallyMap.size();
}

string Float2String(const float flt)
{
    std::ostringstream buff;
    buff << flt;
    return buff.str();
}


///  PERFORM OPERATION
/// //////////////////////
void MeshTallyReader::PerformOperation(string strTally1, string operation, string strTally2)
{
    if(strTally2 == "0")
    {
        cout << "Invalid Input: make sure to pass THREE arguments of types: float char float\n";
        return;
    }

    boost::algorithm::trim(operation);
    boost::algorithm::trim(strTally1);
    boost::algorithm::trim(strTally2);

    // extract integer values from nbTallies
    int talNb1(0), talNb2(0), energyNb1(0), energyNb2(0);

    boost::char_separator<char> sep(".");
    tokenizer tok1(strTally1, sep);
    tokenizer tok2(strTally2, sep);

    //cout << "tok1 -- tok2 " << strTally1 << " -- " << strTally2 << endl; exit(0);

    tokenizer::iterator it = tok1.begin();

    talNb1 = atoi((*it).c_str());
    if(it != tok1.end())
    {
        it++;
        energyNb1 = atoi((*it).c_str());
    }
    it = tok2.begin();
    talNb2 = atoi((*it).c_str());
    if(it != tok2.end())
    {
        it++;
        energyNb2 = atoi((*it).c_str());
    }

    McnpTally tal1, tal2;
    McnpTally newTal;

    for(map<int, McnpTally>::const_iterator it = myTallyMap.begin(); it != myTallyMap.end(); it++)
    {
        McnpTally curTal = (*it).second;

       // cout << "TallyNumbers : " << talNb1 << " -- " << talNb2 << endl; exit(0);

        if(curTal.GetTallyNumber() == talNb1)
            tal1 = curTal;
        if(curTal.GetTallyNumber() == talNb2)
            tal2 = curTal;
    }

    // init new McnpTally
    newTal.AddEnergyTick(0.f);
    newTal.AddEnergyTick(1.f);

    for(list<float>::const_iterator it=tal1.GetXBinBegin(); it!=tal1.GetXBinEnd(); it++)
        newTal.AddXBinTick(*it);
    for(list<float>::const_iterator it=tal1.GetYBinBegin(); it!=tal1.GetYBinEnd(); it++)
        newTal.AddYBinTick(*it);
    for(list<float>::const_iterator it=tal1.GetZBinBegin(); it!=tal1.GetZBinEnd(); it++)
        newTal.AddZBinTick(*it);

    newTal.InitValueArrays();

    // perform operation
    for(int k=0; k<tal1.ZDimension(); k++)
    {
        for(int j=0; j<tal1.YDimension(); j++)
        {
            for(int i=0; i<tal1.XDimension(); i++)
            {
                float result(0.0f);
                float error(0.0f);
                float x, dx, y, dy; // Function G = x{+,-,*,/}y -- with relative errors dx, dy

                x  = tal1.GetValue(energyNb1, i, j, k);
                dx = tal1.GetError(energyNb1, i, j, k);
                y  = tal2.GetValue(energyNb2, i, j, k);
                dy = tal2.GetError(energyNb2, i, j, k);

                if(operation == "+")
                {
                    result = x + y;

                    // Relative Error
                    if(result == 0.f)
                        error = 0.f;
                    else
                        error = (sqrt( x*dx * x*dx + y*dy * y*dy )) / result;
                }
                else if(operation == "-")
                {
                    result = x - y;

                    // Relative Error
                    if(result == 0.f)
                        error = 0.f;
                    else
                        error = (sqrt( x*dx * x*dx + y*dy * y*dy )) / result;

                }
                else if(operation == "/")
                {
                    if(y == 0.0f)
                    {
                        result = 0.0f;
                        error = 0.f;
                    }
                    else
                    {
                        result =  x / y;
                        if(result == 0.f)
                            error = 0.f;
                        else
                        {
                            if(y==0)
                                error = 0.f;
                            else
                                error = sqrt( (dx*x *dx*x) / y / y  + ( x * y*dy * x * y*dy ) / y / y / y / y ) / result;

                        }
                    }
                }
                else if(operation == "*")
                {
                    result = x * y;
                    if(result == 0.f)
                        error = 0.f;
                    else
                        error = sqrt(( x * y*dy * x * y*dy ) +
                                    ( y * x*dx * y * x*dx ));
                }
                newTal.SetValue(0, i, j, k, result);
                newTal.SetError(0, i, j, k, error);
            }
        }
    }

    // find unused tallynumber
    int newTalNb(99);

    for(map<int, McnpTally>::const_iterator it = myTallyMap.begin(); it != myTallyMap.end(); it++)
    {
       McnpTally curTal = (*it).second;
       if(curTal.GetTallyNumber() > newTalNb)
           newTalNb = curTal.GetTallyNumber();
    }

    newTalNb++;
    cout << "New Tally Number : " << newTalNb << endl;
    newTal.SetTallyNumber(newTalNb);

    myTallyMap[newTalNb] = newTal;
}



/// PRINT TALLY INFORMATION
/// /////////////////////////
void MeshTallyReader::PrintTallyInformation() const
{
    string logFileName = myInFileName;
    logFileName.append(".log");
    ofstream log(logFileName.c_str());

    cout << "\n\nTALLY DATA\n=============================\n";
    cout << "Number of Tallies : " << myTallyMap.size() << endl;
    cout << "=============================\n\n";
    log << "\nTALLY DATA\n=============================\n";
    log << "Number of Tallies : " << myTallyMap.size() << endl;
    log << "=============================\n\n";
    for(map<int, McnpTally>::const_iterator it = myTallyMap.begin(); it != myTallyMap.end(); it++)
    {
        McnpTally tally = (*it).second;

        int xSize(1), ySize(1), zSize(1);
        xSize = tally.XDimension()-1;
        ySize = tally.YDimension()-1;
        zSize = tally.ZDimension()-1;

        if(xSize < 1)
            xSize = 1;
        if(ySize < 1)
            ySize = 1;
        if(zSize < 1)
            zSize = 1;

        cout << "Tally Number : " << tally.GetTallyNumber() << endl;
        cout << "---------------------\n";
        cout << "Coordinate System :\n";
        if(tally.GetMeshGeometry() == Mcnp_CARTESIAN)
            cout << "  Cartesian\n";
        else
            cout << "  Cylindrical\n";

        cout << "Energy Bins : \n";
        cout << "  [";
        int numEnergyBins = tally.EnergyDimension();
        /*if(numEnergyBins > 2)
            numEnergyBins--;*/
        for(int e=0; e<numEnergyBins; e++)
        {
            float energyTick = tally.GetEnergyTick(e);
            if(energyTick < 0)
                cout << "Total";
            else
                cout << energyTick;

            if(e == numEnergyBins-1)
                break;

            cout << ", ";
        }
        cout << "]\n";
        cout << "Particle Type : \n";
        if(tally.GetParticleType() == Mcnp_NEUTRON)
            cout << "  Neutron\n";
        else if(tally.GetParticleType() == Mcnp_PHOTON)
            cout << "  Photon\n";
        else
            cout << "unidentified particle type\n";

        cout << "Dimensions of mesh (x,y,z) : \n";
        cout << "  " << tally.XDimension() << " " <<
                        tally.YDimension() << " " <<
                        tally.ZDimension() << endl;

        cout << "Mesh Format : \n";
        if(tally.GetMeshFormat() == Mcnp_IJ)
            cout << "  IJ\n";
        else if(tally.GetMeshFormat() == Mcnp_IK)
            cout << "  IK\n";
        else if(tally.GetMeshFormat() == Mcnp_JK)
            cout << "  JK\n";
        else
            cout << "  COL\n";
        cout <<  "\n\n" << flush;


        log << "Tally Number : " << tally.GetTallyNumber() << endl;
        log << "---------------------\n";
        log << "Coordinate System :\n";
        if(tally.GetMeshGeometry() == Mcnp_CARTESIAN)
            log << "  Cartesian\n";
        else
            log << "  Cylindrical\n";

        log << "Energy Bins : \n";
        log << "  [";
        for(int e=0; e<numEnergyBins; e++)
        {
            float energyTick = tally.GetEnergyTick(e);
            if(energyTick < 0)
                log << "Total";
            else
                log << energyTick;

            if(e == numEnergyBins-1)
                break;

            log << ", ";
        }
        log << "]\n";

        log << "Particle Type : \n";
        if(tally.GetParticleType() == Mcnp_NEUTRON)
            log << "  Neutron\n";
        else if(tally.GetParticleType() == Mcnp_PHOTON)
            log << "  Photon\n";
        else
            log << "unidentified particle type\n";

        log << "Dimensions of mesh (x,y,z) : \n";
        log << "  " <<  tally.XDimension() << " " <<
                        tally.YDimension() << " " <<
                        tally.ZDimension() << endl;

        log << "Mesh Format : \n";
        if(tally.GetMeshFormat() == Mcnp_IJ)
            log << "  IJ\n";
        else if(tally.GetMeshFormat() == Mcnp_IK)
            log << "  IK\n";
        else if(tally.GetMeshFormat() == Mcnp_JK)
            log << "  JK\n";
        else
            log << "  COL\n";
        log <<  "\n\n" << flush;
    }
}


/// PRINT MESHES AVAILABLE
/// ///////////////////////
void MeshTallyReader::PrintMeshesAvailable() const
{
    for(map<int,McnpTally>::const_iterator it = myTallyMap.begin(); it!=myTallyMap.end(); it++)
    {
        McnpTally tally = (*it).second;

        if(tally.EnergyDimension()-1 == 1)
            cout << tally.GetTallyNumber() << "  ";
        else
        {
            for(int i=0; i<tally.EnergyDimension()-1; i++)
                cout << tally.GetTallyNumber() << "." << i << "  ";
        }
        // cout << endl;
    }

    cout << endl;
}


/// PRINT TALLIES AVAILABLE
/// ////////////////////////
void MeshTallyReader::PrintTalliesAvailable()
{
    for(map<int,McnpTally>::const_iterator it = myTallyMap.begin(); it!=myTallyMap.end(); it++)
    {
        McnpTally tally = (*it).second;
        cout << tally.GetTallyNumber() << " ";
    }
    cout << endl;
}



///  READ TALLY FILE
/// //////////////////////
bool MeshTallyReader::ReadTallyFile(string& fileName)
{
    myInFileName = fileName;
    return ReadTallyFile();
}

///  READ TALLY FILE
/// //////////////////////
bool MeshTallyReader::ReadTallyFile()
{
    if(myInFileName == "_initial_")
    {
        cout << "_#_MeshTallyReader::ReadTallyFile() -- Error! No file name specified\n\n ";
        return false;
    }

    cout << "\n--------------------------------------------------------\n"
         << "Reading file : " << myInFileName << "\n--------------------------------------------------------\n" << endl;

  // Open File for reading
  ///////////////////////////
    ifstream readFile;
    readFile.open(myInFileName.c_str(), ios::in);

	if(!readFile.is_open())
	{
        cout << "File " << myInFileName.c_str() << " could not be read!!!\n";
        return false;
	}

  // Read File
  ///////////////////////////////////////
    while(!readFile.eof())
    {
      // get tally number
        int newTallyNumber = GetTallyNumber(readFile);

        if(newTallyNumber == -1) // returns -1 if no new tally number was found, i.e. eof
            break;

        McnpTally currentMcnpTally; // new tally data set
        currentMcnpTally.SetMeshGeometry(Mcnp_CARTESIAN);
        currentMcnpTally.SetTallyNumber(newTallyNumber);

        cout << "\ntally ... " << newTallyNumber << endl;

     // get particle type
        currentMcnpTally.SetParticleType(GetParticleType(readFile));

     // get bins
        GetBinBoundaries(readFile, currentMcnpTally);
        cout << "Dimensions : " << currentMcnpTally.XDimension() << " ... " << currentMcnpTally.YDimension() << " ... " << currentMcnpTally.ZDimension() << endl;

     // get energy bin boundaries
        GetEnergyBinBoundaries(readFile, currentMcnpTally);
        cout << "Number of Energy Bins : " << currentMcnpTally.EnergyDimension() - 1 << endl;

     // determine tally type and read Tally Values
        if(DetermineTallyFormat(readFile, currentMcnpTally) == Mcnp_COL)
            GetColoumnData(readFile, currentMcnpTally);
        else // IJ, IK, JK
            GetMatrixData(readFile, currentMcnpTally);

     // bind tally to tally map
        myTallyMap[currentMcnpTally.GetTallyNumber()] = currentMcnpTally;
    }

    /*if(!CheckTallies())
        cout << "READ FUNCTION NOT VALID\n\n";*/

    PrintTallyInformation();
  //  MakeIntegralList();

    return true;
}

/// SCALE GEOMETRY BY FACTOR
//////////////////////////////
void MeshTallyReader::SetGeometryScalingFactor(const float &scalingFactor)
{
    myGeometryScalingFactor = scalingFactor;
}

/// GET GEOMETRY SCALING FACTOR
/////////////////////////////////
float MeshTallyReader::GetGeometryScalingFactor() const
{
    return myGeometryScalingFactor;
}



/// SET OUTFILE NAME
/////////////////////////////
void MeshTallyReader::SetOutFileName(string &fileName)
{
    if(fileName.find(".vtk") > 0 && fileName.find(".vtk") < fileName.size())
        boost::replace_last(fileName, ".vtk", "_");
    else
        fileName.append("_");

    myOutFileName = fileName;
}


/// WRITE ALL
/////////////////////////////
bool MeshTallyReader::WriteAll()
{
    for(map<int, McnpTally>::iterator it = myTallyMap.begin(); it != myTallyMap.end(); it++) // Tallies
    {
        if(!Write((*it).second.GetTallyNumber(), Write_DATA))
            return false;
        if(!Write((*it).second.GetTallyNumber(), Write_ERROR))
            return false;
    }
    return true;
}


/// WRITE VTK
// print all meshes to files outFileName1.vtk outFileName2.vtk ... outFileNameN.vtk
////////////////////////////////////////////////////////////////////////////////////
bool MeshTallyReader::WriteVTK()
{
    for(map<int, McnpTally>::iterator it = myTallyMap.begin(); it != myTallyMap.end(); it++) // Tallies
        if(!Write((*it).second.GetTallyNumber(), Write_DATA))
            return false;

    return true;
}

/// WRITE VTK
///////////////////////////
bool MeshTallyReader::WriteVTK( string& outName )
{
    myOutFileName = outName;
    return WriteVTK();
}


/// Write Error VTK
////////////////////////////
bool MeshTallyReader::WriteErrorVTK()
{
    for(map<int, McnpTally>::iterator it = myTallyMap.begin(); it != myTallyMap.end(); it++) // Tallies
        if(!Write((*it).second.GetTallyNumber(), Write_ERROR))
            return false;
    return true;
}

/// Write Error VTK
////////////////////////////
bool MeshTallyReader::WriteErrorVTK( string& outName )
{
    myOutFileName = outName;
    return WriteErrorVTK();
}


// Write Selected Tallies
////////////////////////////
bool MeshTallyReader::WriteSelectedTallies(string& selectionList)
{
    boost::char_separator<char> sep(" ,\t");
    tokenizer tok(selectionList, sep);

    for(tokenizer::iterator it = tok.begin(); it != tok.end(); it++)
    {
        // string strTallyNumber(*it);
        int tallyNumber;
        stringstream(*it) >> tallyNumber;
        Write(tallyNumber, Write_DATA);
        Write(tallyNumber, Write_ERROR);
    }

    return true;
}


 /// ***************************************************************
 /// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
 ///                           P R I V A T E
 /// . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
 /// ***************************************************************

/// DETERMINE TALLY FORMAT
///////////////////////////////
McnpMeshFormat MeshTallyReader::DetermineTallyFormat(ifstream &fileStream, McnpTally& curTal)
 {
     string line = GetLine(fileStream);

     if(line.find("Result") != string::npos)
     {
         if(line.find("Energy") != string::npos)
            curTal.SetMeshFormat(Mcnp_COL_WITH_ENERGY);
         else
            curTal.SetMeshFormat(Mcnp_COL);
         return Mcnp_COL;
     }
     else
     {
         // if multiple energy bins -> find != string::npos, we have don't have to skip that line
         if(line.find(" bin") == string::npos)
            line = GetLine(fileStream);

         if(( line.find("Z bin") != string::npos ) ||
            ( line.find("T bin") != string::npos ) )
         {
             if(line.find("T bin") != string::npos )
                 curTal.SetMeshGeometry(Mcnp_CYLINDICAL);
             curTal.SetMeshFormat(Mcnp_IJ);
             return Mcnp_IJ;
         }
         else if(( line.find("Y bin") != string::npos ) ||
                 ( line.find("Z bin") != string::npos ) )
         {
             if(line.find("Z bin") != string::npos )
                 curTal.SetMeshGeometry(Mcnp_CYLINDICAL);
             curTal.SetMeshFormat(Mcnp_IK);
             return Mcnp_IK;
         }
         else
         {
             if(line.find("R bin") != string::npos )
                 curTal.SetMeshGeometry(Mcnp_CYLINDICAL);
             curTal.SetMeshFormat(Mcnp_JK);
             return Mcnp_JK;
         }
     }

     return Mcnp_COL;
 }

/// GET BIN BOUNDARIES
//////////////////////////////
void MeshTallyReader::GetBinBoundaries(ifstream &fileStream, McnpTally& curTalDat)
{
    string line;
    while(!fileStream.eof())
    {
        line = GetLine(fileStream);
        if(line.find("bin boundaries:") != string::npos)
            break;
    }

    /*
      put bin centers in bin lists
    */

    // x bins
    line = GetLine(fileStream);
    {
        bool first(true);
        float previous(0.0);

        boost::char_separator<char> sep("\t ");
        tokenizer tok(line, sep);
        for(tokenizer::iterator it = tok.begin(); it != tok.end(); it++)
        {
            float value(0.0);
            string cur = *it;

            if(!StringToFloat(cur, value))
                continue;

            if(first)
            {
                first=false;
                previous = value;
                continue;
            }

            curTalDat.AddXBinTick((previous+value)/2.0f);
            previous = value;
        }
    }

    // y bins
    line = GetLine(fileStream);
    {
        bool first(true);
        float previous(0.0);

        boost::char_separator<char> sep("\t ");
        tokenizer tok(line, sep);
        for(tokenizer::iterator it = tok.begin(); it != tok.end(); it++)
        {
            float value(0.0);
            string cur = *it;

            if(!StringToFloat(cur, value))
                continue;

            if(first)
            {
                first=false;
                previous = value;
                continue;
            }

            curTalDat.AddYBinTick((previous+value)/2.0f);
            previous = value;
        }
    }

    // z bins
    line = GetLine(fileStream);
    {
        bool first(true);
        float previous(0.0);

        boost::char_separator<char> sep("\t ");
        tokenizer tok(line, sep);
        for(tokenizer::iterator it = tok.begin(); it != tok.end(); it++)
        {
            float value(0.0);
            string cur = *it;

            if(!StringToFloat(cur, value))
                continue;

            if(first)
            {
                first=false;
                previous = value;
                continue;
            }

            curTalDat.AddZBinTick((previous+value)/2.0f);
            previous = value;
        }
    }
}

/// GET COLOUMN DATA
/////////////////////////
void MeshTallyReader::GetColoumnData(ifstream &fileStream, McnpTally &curTal)
{
    string line;
    boost::char_separator<char> sep("\t ");

    // determine which output format we have
    // 1. energy is listed  : Energy X Y Z Value Error
    // 2. energy not listed : X Y Z Value Error
    // we have to check wether the total energy contribution is added after
    // all energy ranges have been covered

    bool haveEnergy(false);

    if(curTal.GetMeshFormat() == Mcnp_COL_WITH_ENERGY)
        haveEnergy = true;

	int progress_N = ( curTal.EnergyDimension() - 1 ) * curTal.XDimension();

    for(int e=0; e<curTal.EnergyDimension()-1; e++)
    {
        float energy = curTal.GetEnergyTick(e+1);

       cout << "\nENERGY : ";
       if(energy < 0 )
           cout << "Total\n";
       else
           cout << energy << endl;

        for(int i=0; i<curTal.XDimension(); i++)
        {
			int progress_n = e * curTal.XDimension() + i + 1;
			PrintProgress(progress_n, progress_N);
            
			for(int j=0; j<curTal.YDimension(); j++)
            {
                for(int k=0; k<curTal.ZDimension(); k++)
                {
                    line = GetLine(fileStream);

                    tokenizer tok(line, sep);

                    tokenizer::iterator it = tok.begin(); // Energy OR x

                    if(haveEnergy || (*it).find("Total") != string::npos)
                        it++; // X
                    it++; // Y
                    it++; // Z
                    it++; // result
                    string resultString = *it;

                    float result;
                    StringToFloat(resultString, result);
                    curTal.SetValue(e, i, j, k, result);

                    if(result < curTal.GetMinVal(energy))
                        curTal.SetMinVal(energy, result);
                    if(result > curTal.GetMaxVal(energy))
                        curTal.SetMaxVal(energy, result);

                    it++;

                    if(it != tok.end()) // errors
                    {
                        string errorString = *it;

                        float error;
                        StringToFloat(errorString, error);
                        curTal.SetError(e, i, j, k, error);
                    }
                }
            }
        }
    }
}

/// GET ENERGY BOUNDARIES
/////////////////////////////
void MeshTallyReader::GetEnergyBinBoundaries(ifstream &fileStream, McnpTally &curTalDat)
{ 
    string line;

    while(!fileStream.eof())
    {
       line = GetLine(fileStream);
       if( line.find("Energy bin boundaries") != string::npos )
           break;
    }

    boost::char_separator<char> sep("\t ");
    tokenizer tok(line, sep);

    float val(0.0);

    for(tokenizer::iterator it = tok.begin(); it != tok.end(); it++)
    {
        string cur = *it;

        if(!StringToFloat(cur, val))
            continue;

        curTalDat.AddEnergyTick(val);
    }

    // now we have all neccessary information available to initialize the dynamic array
    curTalDat.InitValueArrays();
}

/// GET LINE
////////////////////////////
string MeshTallyReader::GetLine(ifstream& fileStream)
{
     char lineIn[65536];
     string line;

     while(!fileStream.eof())
     {
         fileStream.getline(lineIn, 65536);
         line = string(lineIn);
         myLineCount++;

         boost::algorithm::trim(line);

         if(line.length() > 0)
             break;
     }

     return line;
}

/// GET MATRIX DATA
///////////////////////////
void MeshTallyReader::GetMatrixData(ifstream &fileStream, McnpTally &curTal)
{
     string line;

     // Handle different formats !!!
     ///////////////////////////////////////

     bool isInErrorSection(false);
     bool hasError(true);
     bool firstInError(true);

     if(curTal.GetMeshFormat() == Mcnp_IJ)
     {		 
		 int progress_N = ( curTal.EnergyDimension() - 1 ) * curTal.ZDimension();

         for(int e=0; e<curTal.EnergyDimension()-1; e++) // for all energy bins
         {
             float energy = curTal.GetEnergyTick(e+1);

             for(int k=0; k<curTal.ZDimension(); k++)
             {
                 int progress_n = e * curTal.ZDimension() + k;
				 PrintProgress(progress_n, progress_N);

                 if(!isInErrorSection)
                 {
                     do{
                         line = GetLine(fileStream);
                     } while(line.find("Tally Results:") == string::npos);

                     line = GetLine(fileStream); // skip first  line - bin type and range

                     // value block
                     for(int j=0; j<curTal.YDimension(); j++)
                     {
                         line = GetLine(fileStream);

                         boost::char_separator<char> sep("\t ");
                         tokenizer tok(line, sep);

                         int i=0;
                         bool skipFirst(true);
                         for(tokenizer::iterator it = tok.begin(); it != tok.end(); it++)
                         {
                             if(skipFirst) // first entry is bin value
                             {
                                 skipFirst = false;
                                 continue;
                             }
                             string valStr = *it;

                             boost::algorithm::trim(valStr);
                             float val;
                             StringToFloat(valStr, val);

                             if(val < curTal.GetMinVal(energy))
                                 curTal.SetMinVal(energy, val);
                             if(val > curTal.GetMaxVal(energy))
                                 curTal.SetMaxVal(energy, val);

                             curTal.SetValue(e, i, j, k, val);
                             i++;
                         }
                         if(i != curTal.XDimension())
                             cout << "\nError :: " << curTal.XDimension() << " x values expected, " << i << "were found\n";
                     }

                     // handle errors
                     line = GetLine(fileStream); // line is either 'Relative Errors' or '? bin: ... - ...'; can be skipped in any case
                     if(line.find("Relative Errors") == string::npos )
                         hasError = false;

                     if(hasError)
                     {
                         isInErrorSection = true;
                         k--; // rerun k loop with same value for errors
                     }
                 }
                 else // ERROR SECTION
                 {
                     line = GetLine(fileStream); // skip first line, for it contains the bin values

                     // error block
                     for(int j=0; j<curTal.YDimension(); j++)
                     {
                         line = GetLine(fileStream);

                         boost::char_separator<char> sep("\t ");
                         tokenizer tok(line, sep);

                         int i=0;
                         bool skipFirst(true);
                         for(tokenizer::iterator it = tok.begin(); it != tok.end(); it++)
                         {
                             if(skipFirst) // first entry is bin value
                             {
                                 skipFirst = false;
                                 continue;
                             }
                             string valStr = *it;
                             boost::algorithm::trim(valStr);
                             float error;
                             StringToFloat(valStr, error);

                             curTal.SetError(e, i ,j, k, error);
                             i++;
                         }
                     }

                     isInErrorSection = false;
                 }
             }
        }
     }
     else if (curTal.GetMeshFormat() == Mcnp_IK)
     {
		 int progress_N = ( curTal.EnergyDimension() - 1 ) * curTal.YDimension();

         for(int e=0; e<curTal.EnergyDimension()-1; e++) // for all energy bins
         {
             float energy = curTal.GetEnergyTick(e);

             for(int j=0; j<curTal.YDimension(); j++)
             {
				 int progress_n = (e-1) * curTal.ZDimension() + j;
				 PrintProgress(progress_n, progress_N);

                 if(!isInErrorSection)
                 {
                     do{
                         line = GetLine(fileStream);
                     } while(line.find("Tally Results:") == string::npos);

                     line = GetLine(fileStream); // skip first  line - bin type and range

                     // value block
                     for(int k=0; k<curTal.ZDimension(); k++)
                     {
                         line = GetLine(fileStream);

                         boost::char_separator<char> sep("\t ");
                         tokenizer tok(line, sep);

                         int i=0;
                         bool skipFirst(true);
                         for(tokenizer::iterator it = tok.begin(); it != tok.end(); it++)
                         {
                             if(skipFirst) // first entry is bin value
                             {
                                 skipFirst = false;
                                 continue;
                             }
                             string valStr = *it;
                             boost::algorithm::trim(valStr);
                             float val;
                             StringToFloat(valStr, val);

                             if(val < curTal.GetMinVal(energy))
                                 curTal.SetMinVal(energy, val);
                             if(val > curTal.GetMaxVal(energy))
                                 curTal.SetMaxVal(energy, val);

                             curTal.SetValue(e, i, j, k, val);
                             i++;
                         }
                     }

                     // handle errors
                     line = GetLine(fileStream); // line is either 'Relative Errors' or '? bin: ... - ...'; can be skipped in any case
                     if(line.find("Relative Errors") == string::npos )
                         hasError = false;

                     if(hasError)
                     {
                         isInErrorSection = true;
                         j--; // rerun j loop with same value for errors
                     }
                 }
                 else // ERROR SECTION
                 {
                     line = GetLine(fileStream); // skip first line, for it contains the bin values

                     if(firstInError)
                         firstInError = false;
                     else
                         line = GetLine(fileStream); // take second line

                     // error block
                     for(int k=0; k<curTal.ZDimension(); k++)
                     {
                         line = GetLine(fileStream);

                         boost::char_separator<char> sep("\t ");
                         tokenizer tok(line, sep);

                         int i=0;
                         bool skipFirst(true);
                         for(tokenizer::iterator it = tok.begin(); it != tok.end(); it++)
                         {
                             if(skipFirst) // first entry is bin value
                             {
                                 skipFirst = false;
                                 continue;
                             }
                             string valStr = *it;
                             boost::algorithm::trim(valStr);
                             float error;
                             StringToFloat(valStr, error);
                             curTal.SetError(e, i, j, k, error);
                             i++;
                         }
                     }

                     isInErrorSection = false;
                 }
             }
         }
     }
     else // JK
     {
		 int progress_N = ( curTal.EnergyDimension() - 1 ) * curTal.XDimension();
         
         for(int e=0; e<curTal.EnergyDimension()-1; e++) // for all energy bins
         {
             float energy = curTal.GetEnergyTick(e);

             for(int i=0; i<curTal.XDimension(); i++)
             {
				 int progress_n = (e-1) * curTal.ZDimension() + i;
				 PrintProgress(progress_n, progress_N);

                 if(!isInErrorSection)
                 {
                     do{
                         line = GetLine(fileStream);
                     } while(line.find("Tally Results:") == string::npos);

                     line = GetLine(fileStream); // skip second line - bin values

                     // value block
                     for(int k=0; k<curTal.ZDimension(); k++)
                     {
                         line = GetLine(fileStream);

                         boost::char_separator<char> sep("\t ");
                         tokenizer tok(line, sep);

                         int j=0;
                         bool skipFirst(true);
                         for(tokenizer::iterator it = tok.begin(); it != tok.end(); it++)
                         {
                             if(skipFirst) // first entry is bin value
                             {
                                 skipFirst = false;
                                 continue;
                             }
                             string valStr = *it;
                             boost::algorithm::trim(valStr);
                             float val;
                             StringToFloat(valStr, val);

                             if(val < curTal.GetMinVal(energy))
                                 curTal.SetMinVal(energy, val);
                             if(val > curTal.GetMaxVal(energy))
                                 curTal.SetMaxVal(energy, val);

                             curTal.SetValue(e, i, j, k, val);
                             j++;
                         }
                     }

                     // handle errors
                     // if(i==0)
                     {
                         line = GetLine(fileStream); // line is whether 'Relative Errors' or '? bin: ... - ...'; can be skipped in any case

                         if(line.find("Relative Errors") == string::npos )
                             hasError = false;
                     }

                     if(hasError)
                     {
                         isInErrorSection = true;
                         i--; // rerun i loop with same value for errors
                     }
                 }
                 else // ERROR SECTION
                 {
                     line = GetLine(fileStream); // skip first line, for it contains the bin values

                     if(firstInError)
                         firstInError = false;
                     else
                         line = GetLine(fileStream); // take second line

                     // error block
                     for(int k=0; k<curTal.ZDimension(); k++)
                     {
                         line = GetLine(fileStream);

                         boost::char_separator<char> sep("\t ");
                         tokenizer tok(line, sep);

                         int j=0;
                         bool skipFirst(true);
                         for(tokenizer::iterator it = tok.begin(); it != tok.end(); it++)
                         {
                             if(skipFirst) // first entry is bin value
                             {
                                 skipFirst = false;
                                 continue;
                             }
                             string valStr = *it;
                             boost::algorithm::trim(valStr);
                             float error;
                             StringToFloat(valStr, error);
                             curTal.SetError(e, i, j, k, error);
                             j++;
                         }
                     }

                     isInErrorSection = false;
                 }
             }
         }
     }
 }

/// GET PARTICLE TYPE
//////////////////////////
McnpParticleType MeshTallyReader::GetParticleType(ifstream &fileStream)
{
    string line;
    while(!fileStream.eof())
    {
        line = GetLine(fileStream);
// qiu adapte to MCNP6 mesh tally
//		if(line.find("This is a") != string::npos)
        if(line.find("This is a") != string::npos || line.find("mesh tally.") != string::npos)
            break;
    }
// qiu adapte to MCNP6 mesh tally	
/*
    boost::tokenizer<> tok(line);

    string strParticleType;

    int i=1;
    for( boost::tokenizer<>::iterator it = tok.begin(); it!=tok.end(); it++)
    {
        i++;
        strParticleType = *it;

        if(i==5)
            break;
    }
*/
    McnpParticleType particleType;
// qiu adapte to MCNP6 mesh tally
//    if(strParticleType == "neutron")
    if(line.find("neutron") != string::npos)
    {
        particleType = Mcnp_NEUTRON;
        cout << "Neutron Tally\n";
    }
//    else if(strParticleType == "photon")
    else if(line.find("photon") != string::npos)
    {
        particleType = Mcnp_PHOTON;
        cout << "Photon Tally\n";
    }
    else
        cout << "Undefined Particle Type\n";

    return particleType;
}


/// GET TALLY
/// /////////////////////
McnpTally MeshTallyReader::GetTally(int tallyNumber)
{
    McnpTally returnTally = myTallyMap[tallyNumber];
    return returnTally;
}

/// GET TALLY NUMBER
/////////////////////////
int MeshTallyReader::GetTallyNumber(ifstream &fileStream)
{
    int tallyNumber(-1);

    while(true && !fileStream.eof())
    {
        string line = GetLine(fileStream);

        if(line.find("Mesh Tally Number") == string::npos ) // skip all lines that do not contain "Mesh Tally Number"
            continue;

        boost::tokenizer<> tok(line);
        string last;

        for( boost::tokenizer<>::iterator it = tok.begin(); it!=tok.end(); it++)
           last = *it;

        tallyNumber = boost::lexical_cast<int> (last);
        break;
    }

    return tallyNumber;
}

///   INIT
////////////////////////////
void MeshTallyReader::Init()
 {
    myGeometryScalingFactor = 1.0;
    myInFileName = "_initial_";
 }

/// STRING TO float
///////////////////////////
bool MeshTallyReader::StringToFloat(string str, float& val)
{
    stringstream ss(str);

    if( (ss >> val).fail() )
        return false;

    return true;
}


/// PRINT PROGRESS
///////////////////////////
void MeshTallyReader::PrintProgress(const int& n, const int& N) const
{
	float progress = float(n) / float(N) * 100.f;
	cout.precision(2);
	cout.setf(ios::fixed);
        cout << setw(5) << "\r" << progress << "%" << flush;
        cout.unsetf(ios::fixed);
}


/// WRITE
///////////////////////////
bool MeshTallyReader::Write(int tallyNumber, WriteType wType)
{
  //  for(map<int, McnpTally>::iterator it = myTallyMap.begin(); it != myTallyMap.end(); it++) // Tallies
    {
        McnpTally curTal = myTallyMap[tallyNumber]; //(*it).second;

        int dimX = curTal.XDimension();
        int dimY = curTal.YDimension();
        int dimZ = curTal.ZDimension();

        int energyBins = curTal.EnergyDimension()-1;

        for(int e=0; e<energyBins; e++) // Energies
        {

            // set output name
            string outFileName = myOutFileName;

            stringstream floatToString;
            floatToString << curTal.GetTallyNumber();
            outFileName.append(floatToString.str());

            if(curTal.EnergyDimension() > 2) // multiple energy bins - output format: FILENAME_TALLYNUMBER.ENERGYNUMBER
            {
                stringstream energyToString;
                energyToString << e;
                outFileName.append(".");
                outFileName.append(energyToString.str());
            }

            if(wType == Write_ERROR)
                outFileName.append("_error");

            outFileName.append(".vtk");

            cout << " writing ";
            if(wType == Write_ERROR)
                cout << "error ";
            cout <<	"vtk file : " << outFileName << endl;

            ofstream outFile(outFileName.c_str());

            // write mesh
	    ////////////////////
            outFile << "# vtk DataFile Version 3.0\n"
                    << "Karlsruhe Institute of Technology - mt2vtk\n"
                    << "ASCII\n"
                    << "DATASET STRUCTURED_GRID\n"
                    << "DIMENSIONS " << dimX << " " << dimY << " " << dimZ << endl
                    << "POINTS " << dimX * dimY * dimZ << " float\n";

            int lineBreakCount(0);
            for(list<float>::const_iterator itZ=curTal.GetZBinBegin(); itZ!=curTal.GetZBinEnd(); itZ++)
            {
                for(list<float>::const_iterator itY=curTal.GetYBinBegin(); itY!=curTal.GetYBinEnd(); itY++)
                {
                    for(list<float>::const_iterator itX=curTal.GetXBinBegin(); itX!=curTal.GetXBinEnd(); itX++)
                    {
                        if(++lineBreakCount > 3)
                        {
                            lineBreakCount = 1;
                            outFile << endl;
                        }
                        outFile << myGeometryScalingFactor * (*itX) << " " <<
                                   myGeometryScalingFactor * (*itY) << " " <<
                                   myGeometryScalingFactor * (*itZ) << " ";
                    }
                }
            }

            // write values
	    /////////////////////
            outFile << "\nPOINT_DATA " << curTal.XDimension() * curTal.YDimension() * curTal.ZDimension() << endl
                    << "SCALARS scalars float\n"
                    << "LOOKUP_TABLE default\n";

            int kOffset, jOffset, offset;
            int i(0), j(0), k(0);

            lineBreakCount = 0;
            for(list<float>::const_iterator itZ=curTal.GetZBinBegin(); itZ!=curTal.GetZBinEnd(); itZ++)
            {
                kOffset = k * dimX * dimY;
                j=0;

                for(list<float>::const_iterator itY=curTal.GetYBinBegin(); itY!=curTal.GetYBinEnd(); itY++)
                {
                    jOffset = j * dimX;

                    i=0;

                    for(list<float>::const_iterator itX=curTal.GetXBinBegin(); itX!=curTal.GetXBinEnd(); itX++)
                    {
                        offset = i + jOffset + kOffset;

                        float value;

                        if(wType == Write_ERROR)
                            value = curTal.GetError(e, i, j, k);
                        else
                            value = curTal.GetValue(e, i, j, k);

                        if(++lineBreakCount > 9)
                        {
                            lineBreakCount = 1;
                            outFile << endl;
                        }

                        outFile << value << " ";

                        i++;
                    }
                    j++;
                }
                k++;
            }
        }
    }

    cout << "... done\n";

    return true;
}
