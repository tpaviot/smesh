#include <math.h>

#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeSphere.hxx>
#include <BRepPrimAPI_MakeTorus.hxx>
#include <SMESH_Gen.hxx>
#include <StdMeshers_AutomaticLength.hxx>
#include <StdMeshers_TrianglePreference.hxx>
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
    StdMeshers_AutomaticLength* hyp1d = new StdMeshers_AutomaticLength(0,0,meshgen);
    StdMeshers_TrianglePreference* hyp2d = new StdMeshers_TrianglePreference(1,0,meshgen);
    mesh->AddHypothesis(mesh->GetShapeToMesh(), 0);
    mesh->AddHypothesis(mesh->GetShapeToMesh(), 1);
    // compute the mesh
    meshgen->Compute(*mesh, mesh->GetShapeToMesh());
    // free memory
    delete meshgen;
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
    delete meshgen;
    delete mesh;
}

TEST(MeshBasicGeometriesSuite, testMeshTorus)
{
    // the same as the previous test, but with a sphere
    BRepPrimAPI_MakeTorus my_torus(10., 20.);
    my_torus.Build();
    ASSERT_TRUE(my_torus.IsDone());
    // create the Mesh
    SMESH_Gen* meshgen = new SMESH_Gen();
    SMESH_Mesh* mesh = meshgen->CreateMesh(0, true);
    // set geometry to be meshed
    mesh->ShapeToMesh(my_torus.Shape());
    ASSERT_TRUE(mesh->HasShapeToMesh());
    // compute the mesh
    meshgen->Compute(*mesh, mesh->GetShapeToMesh());
    // free memory
    delete meshgen;
    delete mesh;
}

int main(int argc, char **argv){
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
