/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkImplicitPolyData.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkImplicitPolyData
// .SECTION Description
//
// Implicit function that computes the distance from a point x to the
// nearest point p on an input vtkPolyData. The sign of the function
// is set to the sign of the dot product between the angle-weighted
// pseudonormal at the nearest surface point and the vector x - p.
// Points interior to the geometry have a negative distance, points on
// the exterior have a positive distance, and points on the input
// vtkPolyData have a distance of zero. The gradient of the function
// is the angle-weighted pseudonormal at the nearest point.
//
// Baerentzen, J. A. and Aanaes, H. (2005). Signed distance
// computation using the angle weighted pseudonormal. IEEE
// Transactions on Visualization and Computer Graphics, 11:243-253.
//
// Written by Chris Weigle and Cory Quammen, The University of North
// Carolina at Chapel Hill.

#ifndef __vtkImplicitPolyData_h
#define __vtkImplicitPolyData_h

#include "vtkImplicitFunction.h"


class vtkCellLocator;
class vtkPolyData;
class vtkTriangleFilter;

class vtkImplicitPolyData : public vtkImplicitFunction
{
public:
  static vtkImplicitPolyData *New();
  vtkTypeMacro(vtkImplicitPolyData,vtkImplicitFunction);
  void PrintSelf(ostream& os, vtkIndent indent);

  // Description:
  // Return the MTime also considering the Input dependency.
  unsigned long GetMTime();

  // Description:
  // Evaluate plane equation of nearest triangle to point x[3].
  double EvaluateFunction(double x[3]);

  // Description:
  // Evaluate function gradient of nearest triangle to point x[3].
  void EvaluateGradient(double x[3], double g[3]);

  // Description:
  // Set the input vtkPolyData used for the implicit function
  // evaluation.  Passes input through an internal instance of
  // vtkTriangleFilter to remove vertices and lines, leaving only
  // triangular polygons for evaluation as implicit planes.
  void SetInput(vtkPolyData *input);

  // Description:
  // Set/get the function value to use if no input vtkPolyData
  // specified.
  vtkSetMacro(NoValue, double);
  vtkGetMacro(NoValue, double);

  // Description:
  // Set/get the function gradient to use if no input vtkPolyData
  // specified.
  vtkSetVector3Macro(NoGradient, double);
  vtkGetVector3Macro(NoGradient, double);

  // Description:
  // Set/get the tolerance usued for the locator.
  vtkGetMacro(Tolerance, double);
  vtkSetMacro(Tolerance, double);

protected:
  vtkImplicitPolyData();
  ~vtkImplicitPolyData();

  double NoValue;
  double NoGradient[3];
  double Tolerance;

  double SharedEvaluate( double x[3], double n[3] );

private:
  vtkImplicitPolyData(const vtkImplicitPolyData&);  // Not implemented.
  void operator=(const vtkImplicitPolyData&);  // Not implemented.

  vtkTriangleFilter *TriangleFilter;
  vtkPolyData       *Input;
  vtkCellLocator    *Locator;

};

#endif
