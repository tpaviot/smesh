#include <SMESH_Gen.hxx>
#include <gtest/gtest.h>

TEST(CreateMeshTestSuite, testCreateSimpleMesh)
{
    SMESH_Gen* meshgen = new SMESH_Gen();
    SMESH_Mesh* mesh = meshgen->CreateMesh(0, true);
    delete meshgen;
}

TEST(CreateMeshTestSuite, testCreate2Meshes)
{
    SMESH_Gen* meshgen = new SMESH_Gen();
    SMESH_Mesh* mesh1 = meshgen->CreateMesh(0, true);
    SMESH_Mesh* mesh2 = meshgen->CreateMesh(1, true);
    delete meshgen;
}

TEST(CreateMeshTestSuite, testCreateSimpleMeshNotEmbeddedMode)
{
    SMESH_Gen* meshgen = new SMESH_Gen();
    SMESH_Mesh* mesh = meshgen->CreateMesh(0, false);
    delete meshgen;
}

TEST(CreateMeshTestSuite, testCreateSimpleMeshWithDefaultSegments)
{
    SMESH_Gen* meshgen = new SMESH_Gen();
    meshgen->SetDefaultNbSegments(300);
    ASSERT_EQ(meshgen->GetDefaultNbSegments(), 300);
}

int main(int argc, char **argv){
  testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}
