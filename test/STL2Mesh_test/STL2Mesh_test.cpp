#include <SMESH_Gen.hxx>
#include <gtest/gtest.h>

#include "STL2Mesh_test_config.h"

TEST(STL2MeshTestSuite, testSTLAscii)
{
    SMESH_Gen* meshgen = new SMESH_Gen();
    SMESH_Mesh* mesh = meshgen->CreateMesh(0, true);
    mesh->STLToMesh(cube_stl_ascii);
    ASSERT_EQ(mesh->NbNodes(), 8);
    ASSERT_EQ(mesh->NbFaces(), 12);
    delete meshgen;
    delete mesh;
}

TEST(STL2MeshTestSuite, testSTLBinary)
{
    SMESH_Gen* meshgen = new SMESH_Gen();
    SMESH_Mesh* mesh = meshgen->CreateMesh(0, true);
    mesh->STLToMesh(cube_stl_binary);
    ASSERT_EQ(mesh->NbNodes(), 8);
    ASSERT_EQ(mesh->NbFaces(), 12);
    delete meshgen;
    delete mesh;
}

int main(int argc, char **argv){
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
