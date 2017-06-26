#include <math.h>

#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepPrimAPI_MakeTorus.hxx>
#include <SMESH_Mesh.hxx>
#include <SMESH_Gen.hxx>
#include <StdMeshers_AutomaticLength.hxx>
#include <StdMeshers_Arithmetic1D.hxx>
#include <StdMeshers_NumberOfSegments.hxx>
#include <StdMeshers_Regular_1D.hxx>
#include <StdMeshers_MEFISTO_2D.hxx>
#include <StdMeshers_Quadrangle_2D.hxx>

#include <gtest/gtest.h>

TEST(MeshBasicGeometriesSuite, testMeshBox)
{
    // create a box, mixing integers and floats
    BRepPrimAPI_MakeBox my_box(10.,10.,10.);
    my_box.Build();
    ASSERT_TRUE(my_box.IsDone());
    // create the Mesh
    SMESH_Gen* meshgen = new SMESH_Gen();
    SMESH_Mesh* mesh = meshgen->CreateMesh(0, true);
    // set geometry to be meshed
    mesh->ShapeToMesh(my_box.Shape());
    ASSERT_TRUE(mesh->HasShapeToMesh());
    // check bounding box. It should be 10.sqrt(3)==17.32050807568877
    double diagonal_size = mesh->GetShapeDiagonalSize(mesh->GetShapeToMesh());
    ASSERT_GT(diagonal_size, 17.320508);
    ASSERT_LT(diagonal_size, 17.320509);
    // create and add hypothesis
    StdMeshers_AutomaticLength* hyp1d_0 = new StdMeshers_AutomaticLength(0,0,meshgen);
    StdMeshers_NumberOfSegments* hyp1d_1 = new StdMeshers_NumberOfSegments(1,0,meshgen);
    hyp1d_1->SetNumberOfSegments(1);
    StdMeshers_Regular_1D* hyp1d_2 = new StdMeshers_Regular_1D(2,0,meshgen);
    mesh->AddHypothesis(mesh->GetShapeToMesh(), 0);
    mesh->AddHypothesis(mesh->GetShapeToMesh(), 1);
    mesh->AddHypothesis(mesh->GetShapeToMesh(), 2);
    // compute the mesh
    meshgen->Compute(*mesh, mesh->GetShapeToMesh());
    // free memory
    delete mesh;
}

TEST(MeshBasicGeometriesSuite, testMeshBoxMEFISTO2)
{
    // create a box, mixing integers and floats
    BRepPrimAPI_MakeBox my_box(10.,10.,10.);
    my_box.Build();
    ASSERT_TRUE(my_box.IsDone());
    // create the Mesh
    SMESH_Gen* meshgen = new SMESH_Gen();
    SMESH_Mesh* mesh = meshgen->CreateMesh(0, true);
    // set geometry to be meshed
    mesh->ShapeToMesh(my_box.Shape());
    ASSERT_TRUE(mesh->HasShapeToMesh());
    // check bounding box. It should be 10.sqrt(3)==17.32050807568877
    double diagonal_size = mesh->GetShapeDiagonalSize(mesh->GetShapeToMesh());
    ASSERT_GT(diagonal_size, 17.320508);
    ASSERT_LT(diagonal_size, 17.320509);
    // create and add hypothesis
    StdMeshers_Arithmetic1D* hyp1d = new StdMeshers_Arithmetic1D(0,0,meshgen);
    hyp1d->SetLength(0.1, false); // the smallest distance between 2 points
    hyp1d->SetLength(0.5, true);  // the longest distance between 2 points
    StdMeshers_Regular_1D* an1DAlgo = new StdMeshers_Regular_1D(1, 0, meshgen); // interpolation
    StdMeshers_MEFISTO_2D* mef2d = new StdMeshers_MEFISTO_2D(2,0,meshgen) ;
    mesh->AddHypothesis(mesh->GetShapeToMesh(), 0);
    mesh->AddHypothesis(mesh->GetShapeToMesh(), 1);
    mesh->AddHypothesis(mesh->GetShapeToMesh(), 2);
    // compute the mesh
    meshgen->Compute(*mesh, mesh->GetShapeToMesh());
    // free memory
    delete mesh;
}

TEST(MeshBasicGeometriesSuite, testMeshSphere)
{
    // the same as the previous test, but with a sphere
    BRepPrimAPI_MakeSphere my_sphere(10.);
    my_sphere.Build();
    ASSERT_TRUE(my_sphere.IsDone());
    // create the Mesh
    SMESH_Gen* meshgen = new SMESH_Gen();
    SMESH_Mesh* mesh = meshgen->CreateMesh(0, true);
    // set geometry to be meshed
    mesh->ShapeToMesh(my_sphere.Shape());
    ASSERT_TRUE(mesh->HasShapeToMesh());
    // compute the mesh
    meshgen->Compute(*mesh, mesh->GetShapeToMesh());
    // free memory
    delete mesh;
}

TEST(MeshBasicGeometriesSuite, testMeshQuadrangleTorus)
{
    // the same as the previous test, but with a sphere
    BRepPrimAPI_MakeTorus my_torus(10., 20.);
    my_torus.Build();
    ASSERT_TRUE(my_torus.IsDone());
    // create the Mesh
    SMESH_Gen* meshgen = new SMESH_Gen();
    SMESH_Mesh* mesh = meshgen->CreateMesh(0, true);
    // 1D
    StdMeshers_Arithmetic1D* hyp1d = new StdMeshers_Arithmetic1D(0,0,meshgen);
    hyp1d->SetLength(0.1, false); // the smallest distance between 2 points
    hyp1d->SetLength(0.5, true);  // the longest distance between 2 points
    StdMeshers_Regular_1D* an1DAlgo = new StdMeshers_Regular_1D(1, 0, meshgen); // interpolation
    StdMeshers_Quadrangle_2D *a2dAlgo = new StdMeshers_Quadrangle_2D(2, 0, meshgen);
    // add hypothesis
    mesh->AddHypothesis(mesh->GetShapeToMesh(), 0);
    mesh->AddHypothesis(mesh->GetShapeToMesh(), 1);
    mesh->AddHypothesis(mesh->GetShapeToMesh(), 2);
    // set geometry to be meshed
    mesh->ShapeToMesh(my_torus.Shape());
    ASSERT_TRUE(mesh->HasShapeToMesh());
    // compute the mesh
    meshgen->Compute(*mesh, mesh->GetShapeToMesh());
    // free memory
    delete mesh;
}

int main(int argc, char **argv){
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
