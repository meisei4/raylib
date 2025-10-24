/*******************************************************************************************
*
*  WIP  raylib [shapes] fixed function didactic example - World -> NDC -> screen
*
*   Example complexity rating: [★★★☆] 3/4 ??????
*
*   Focus: Observed basis concepts, frustum, near-plane intersections, world->NDC remap into
*          isotropic NDC cube, and "intentionally" incorrect near-plane texturing.
*
********************************************************************************************/

#include "raylib.h"
#include "raymath.h"
#include "rlgl.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

//----------------------------------------------------------------------------------
// Constants
//----------------------------------------------------------------------------------
static const int N64_WIDTH  = 640;
static const int N64_HEIGHT = 480;

static const float ANGULAR_VELOCITY = 1.25f;
static const float FOVY_PERSPECTIVE = 60.0f;

static const Vector3 MODEL_POS   = { 0.0f, 0.0f, 0.0f };
static const Vector3 MODEL_SCALE = { 1.0f, 1.0f, 1.0f };

static const Vector3 OBSERVER_POS   = { 0.0f, 0.0f, 2.0f };
static const Vector3 JUGEMU_POS_ISO = { 3.0f, 1.0f, 3.0f };

#define PERSPECTIVE_CORRECT 0
#define DEPTH_TEST_ON 1

//----------------------------------------------------------------------------------
// Small math helpers
//----------------------------------------------------------------------------------
static inline Vector3 observedLineOfSight(Camera3D *observer)
{
    Vector3 lineOfSight = Vector3Subtract(observer->target, observer->position);
    return Vector3Normalize(lineOfSight);
}

static inline void observedBasis(Camera3D *observer, Vector3 *observedLineOfSightOut, Vector3 *observedRightOut, Vector3 *observedUpOut)
{
    Vector3 los = observedLineOfSight(observer);
    Vector3 right = Vector3Normalize(Vector3CrossProduct(los, observer->up));
    Vector3 up = Vector3Normalize(Vector3CrossProduct(right, los));
    *observedLineOfSightOut = los;
    *observedRightOut = right;
    *observedUpOut = up;
}

static inline Vector3 rotatePointAboutAxis(Vector3 point, Vector3 axisStart, Vector3 axisEnd, float angleRadians)
{
    Vector3 axisDir = Vector3Normalize(Vector3Subtract(axisEnd, axisStart));
    Vector3 localFromAxis = Vector3Subtract(point, axisStart);
    Vector3 rotatedLocal = Vector3RotateByAxisAngle(localFromAxis, axisDir, angleRadians);
    return Vector3Add(axisStart, rotatedLocal);
}

static inline void rotateVerticesInPlaneSlice(Vector3 *vertices, int vertexCount, float rotationRadians)
{
    float s = sinf(rotationRadians), c = cosf(rotationRadians);
    for (int vertexIndex = 0; vertexIndex < vertexCount; vertexIndex++)
    {
        float x0 = vertices[vertexIndex].x, z0 = vertices[vertexIndex].z;
        vertices[vertexIndex].x =  c * x0 + s * z0;
        vertices[vertexIndex].z = -s * x0 + c * z0;
    }
}

static inline Vector3 applyModelTranslateRotateScale(Vector3 modelCoord, Vector3 modelPosition, Vector3 modelScale, float rotationRadians)
{
    Matrix M = MatrixMultiply(
                    MatrixMultiply(
                        MatrixScale(modelScale.x, modelScale.y, modelScale.z),
                        MatrixRotateY(rotationRadians)),
                    MatrixTranslate(modelPosition.x, modelPosition.y, modelPosition.z));
    return Vector3Transform(modelCoord, M);
}

static inline Vector3 applyInverseModelTranslateRotateScale(Vector3 worldCoord, Vector3 modelPosition, Vector3 modelScale, float rotationRadians)
{
    Matrix M = MatrixMultiply(
                    MatrixMultiply(
                        MatrixScale(modelScale.x, modelScale.y, modelScale.z),
                        MatrixRotateY(rotationRadians)),
                    MatrixTranslate(modelPosition.x, modelPosition.y, modelPosition.z));
    Matrix invM = MatrixInvert(M);
    return Vector3Transform(worldCoord, invM);
}

static inline Vector3 nearPlaneIntersection(Camera3D *observer, float nearClipPlane, Vector3 worldCoord)
{
    Vector3 los = observedLineOfSight(observer);
    Vector3 rayToWorld = Vector3Subtract(worldCoord, observer->position);
    float signedDepth = Vector3DotProduct(rayToWorld, los); // positive forward depth
    float t = nearClipPlane / (signedDepth + 1e-9f);
    return Vector3Add(observer->position, Vector3Scale(rayToWorld, t));
}

static inline Vector3 triangleNormal(Vector3 a, Vector3 b, Vector3 c)
{
    return Vector3Normalize(Vector3CrossProduct(Vector3Subtract(b, a), Vector3Subtract(c, a)));
}

//----------------------------------------------------------------------------------
// Drawing helpers
//----------------------------------------------------------------------------------
static void drawObservedAxes(Camera3D *observer)
{
    Vector3 los, right, up; observedBasis(observer, &los, &right, &up);
    rlSetLineWidth(1.5f);
    DrawLine3D(observer->position, Vector3Add(observer->position, right), PURPLE);
    DrawLine3D(observer->position, Vector3Add(observer->position, up), DARKGREEN);
    DrawLine3D(observer->position, Vector3Add(observer->position, los), SKYBLUE);
}

static void drawFrustum(Camera3D *observer, float aspect, float nearClipPlane, float farClipPlane)
{
    Vector3 los, right, up; observedBasis(observer, &los, &right, &up);

    Vector3 centerNear = Vector3Add(observer->position, Vector3Scale(los, nearClipPlane));
    Vector3 centerFar  = Vector3Add(observer->position, Vector3Scale(los,  farClipPlane));

    float halfFovy = DEG2RAD * observer->fovy * 0.5f;
    float halfHeightNear = nearClipPlane * tanf(halfFovy);
    float halfWidthNear  = halfHeightNear * aspect;
    float halfHeightFar  = farClipPlane * tanf(halfFovy);
    float halfWidthFar   = halfHeightFar * aspect;

    Vector3 nearTL = Vector3Add(Vector3Add(centerNear, Vector3Scale(up,  halfHeightNear)), Vector3Scale(right, -halfWidthNear));
    Vector3 nearTR = Vector3Add(Vector3Add(centerNear, Vector3Scale(up,  halfHeightNear)), Vector3Scale(right,  halfWidthNear));
    Vector3 nearBR = Vector3Add(Vector3Add(centerNear, Vector3Scale(up, -halfHeightNear)), Vector3Scale(right,  halfWidthNear));
    Vector3 nearBL = Vector3Add(Vector3Add(centerNear, Vector3Scale(up, -halfHeightNear)), Vector3Scale(right, -halfWidthNear));

    Vector3 farTL = Vector3Add(Vector3Add(centerFar, Vector3Scale(up,  halfHeightFar)), Vector3Scale(right, -halfWidthFar));
    Vector3 farTR = Vector3Add(Vector3Add(centerFar, Vector3Scale(up,  halfHeightFar)), Vector3Scale(right,  halfWidthFar));
    Vector3 farBR = Vector3Add(Vector3Add(centerFar, Vector3Scale(up, -halfHeightFar)), Vector3Scale(right,  halfWidthFar));
    Vector3 farBL = Vector3Add(Vector3Add(centerFar, Vector3Scale(up, -halfHeightFar)), Vector3Scale(right, -halfWidthFar));

    rlSetLineWidth(1.0f);
    DrawLine3D(nearTL, nearTR, SKYBLUE); DrawLine3D(nearTR, nearBR, SKYBLUE);
    DrawLine3D(nearBR, nearBL, SKYBLUE); DrawLine3D(nearBL, nearTL, SKYBLUE);

    DrawLine3D(farTL, farTR, GRAY); DrawLine3D(farTR, farBR, GRAY);
    DrawLine3D(farBR, farBL, GRAY); DrawLine3D(farBL, farTL, GRAY);

    DrawLine3D(nearTL, farTL, DARKBLUE); DrawLine3D(nearTR, farTR, DARKBLUE);
    DrawLine3D(nearBR, farBR, DARKBLUE); DrawLine3D(nearBL, farBL, DARKBLUE);

}

static void drawNdcCube(Camera3D *observer,
                        Vector3 centerNear,
                        float halfWidthNear, float halfHeightNear, float halfDepthNdcCube,
                        Color edgeColor)
{
    Vector3 los, right, up; observedBasis(observer, &los, &right, &up);

    Vector3 corners[8];
    int writeIndex = 0;
    for (int zSign = -1; zSign <= 1; zSign += 2)
    for (int ySign = -1; ySign <= 1; ySign += 2)
    for (int xSign = -1; xSign <= 1; xSign += 2)
    {
        corners[writeIndex++] = Vector3Add(centerNear,
             Vector3Add(Vector3Scale(right, xSign*halfWidthNear),
             Vector3Add(Vector3Scale(up,    ySign*halfHeightNear),
                        Vector3Scale(los, zSign*halfDepthNdcCube))));
    }

    Vector3 nearTL = corners[3], nearTR = corners[2], nearBR = corners[0], nearBL = corners[1];
    Vector3  farTL = corners[7],  farTR = corners[6],  farBR = corners[4],  farBL = corners[5];

    rlSetLineWidth(1.0f);
    DrawLine3D(nearTL, nearTR, edgeColor); DrawLine3D(nearTR, nearBR, edgeColor);
    DrawLine3D(nearBR, nearBL, edgeColor); DrawLine3D(nearBL, nearTL, edgeColor);
    DrawLine3D(farTL,  farTR,  edgeColor); DrawLine3D(farTR,  farBR,  edgeColor);
    DrawLine3D(farBR,  farBL,  edgeColor); DrawLine3D(farBL,  farTL,  edgeColor);
    DrawLine3D(nearTL, farTL,  edgeColor); DrawLine3D(nearTR, farTR,  edgeColor);
    DrawLine3D(nearBR, farBR,  edgeColor); DrawLine3D(nearBL, farBL,  edgeColor);
}

//----------------------------------------------------------------------------------
// World -> NDC mapping
//----------------------------------------------------------------------------------
static Vector3 worldCoordToNdcCoord(float aspect,
                                    Camera3D *observer,
                                    float nearClipPlane, float farClipPlane,
                                    Vector3 worldCoord)
{
    Vector3 los, right, up; observedBasis(observer, &los, &right, &up);

    float halfFovy = DEG2RAD * observer->fovy * 0.5f;
    float halfHeightNear = nearClipPlane * tanf(halfFovy);
    float halfWidthNear  = halfHeightNear * aspect;

    Vector3 intersectionCoord = nearPlaneIntersection(observer, nearClipPlane, worldCoord);
    Vector3 centerNear = Vector3Add(observer->position, Vector3Scale(los, nearClipPlane));
    Vector3 clipPlaneVector = Vector3Subtract(intersectionCoord, centerNear);

    float xNdc = Vector3DotProduct(clipPlaneVector, right) / (halfWidthNear + 1e-9f);
    float yNdc = Vector3DotProduct(clipPlaneVector, up)    / (halfHeightNear + 1e-9f);

    float signedDepth = Vector3DotProduct(Vector3Subtract(worldCoord, observer->position), los);
    float zNdc = ((farClipPlane + nearClipPlane) - (2.0f * farClipPlane * nearClipPlane) / (signedDepth + 1e-9f))
               / (farClipPlane - nearClipPlane);

    return (Vector3){ xNdc, yNdc, zNdc };
}

static Vector3 scaleNdcCoordByNearClipPlane(Camera3D *observer,
                                            Vector3 centerNear,
                                            float halfWidthNear, float halfHeightNear, float halfDepthNdcCube,
                                            Vector3 ndcCoord)
{
    Vector3 los, right, up; observedBasis(observer, &los, &right, &up);
    return Vector3Add(centerNear,
           Vector3Add(Vector3Scale(right, ndcCoord.x * halfWidthNear),
           Vector3Add(Vector3Scale(up,    ndcCoord.y * halfHeightNear),
                      Vector3Scale(los,   ndcCoord.z * halfDepthNdcCube))));
}

static void updateWorldToNdcMappedMesh(Mesh *ndcMesh, Mesh *worldMesh,
                                       int screenWidth, int screenHeight,
                                       Camera3D *observer,
                                       float nearClipPlane, float farClipPlane,
                                       Vector3 modelPosition, Vector3 modelScale, float meshRotationRadians)
{
    float aspect = (float)screenWidth / (float)screenHeight;
    Vector3 los = observedLineOfSight(observer);

    float halfFovy = DEG2RAD * observer->fovy * 0.5f;
    float halfHeightNear = nearClipPlane * tanf(halfFovy);

    // --- ISOTROPIC DIDACTIC VIEW ---
    float halfWidthNear  = halfHeightNear;
    float halfDepthNdc   = halfHeightNear;
    // // --- ANISOTROPIC VIEW ---
    // float halfWidthNear = halfHeightNear * aspect;
    // float halfDepthNdc  = 0.5f * (farClipPlane - nearClipPlane);

    Vector3 centerNear = Vector3Add(observer->position, Vector3Scale(los, nearClipPlane));
    Vector3 ndcCubeCenter  = Vector3Add(centerNear, Vector3Scale(los, halfDepthNdc));

    if (worldMesh->indices && ndcMesh->indices)
    {
        memcpy(ndcMesh->indices, worldMesh->indices, sizeof(unsigned short) * 3 * worldMesh->triangleCount);
        ndcMesh->triangleCount = worldMesh->triangleCount;
    }

    for (int vertexIndex = 0; vertexIndex < worldMesh->vertexCount; vertexIndex++)
    {
        Vector3 objectVertex = { worldMesh->vertices[3*vertexIndex+0], worldMesh->vertices[3*vertexIndex+1], worldMesh->vertices[3*vertexIndex+2] };
        Vector3 worldVertex = applyModelTranslateRotateScale(objectVertex, modelPosition, modelScale, meshRotationRadians);
        Vector3 ndcCoord = worldCoordToNdcCoord(aspect, observer, nearClipPlane, farClipPlane, worldVertex);
        Vector3 scaledNdcCoord = scaleNdcCoordByNearClipPlane(observer, ndcCubeCenter, halfWidthNear, halfHeightNear, halfDepthNdc, ndcCoord);
        Vector3 mappedObjectCoord = applyInverseModelTranslateRotateScale(scaledNdcCoord, modelPosition, modelScale, meshRotationRadians);

        ndcMesh->vertices[3*vertexIndex+0] = mappedObjectCoord.x;
        ndcMesh->vertices[3*vertexIndex+1] = mappedObjectCoord.y;
        ndcMesh->vertices[3*vertexIndex+2] = mappedObjectCoord.z;
    }

    UploadMesh(ndcMesh, false);
}

static void mapFrustumToNdcCube(Camera3D *observer,
                                int screenWidth, int screenHeight,
                                float nearClipPlane, float farClipPlane,
                                Model *worldModel, Model *ndcModel,
                                Vector3 modelPosition, Vector3 modelScale, float meshRotationRadians)
{
    Vector3 los = observedLineOfSight(observer);
    float halfFovy = DEG2RAD * observer->fovy * 0.5f;
    float halfHeightNear = nearClipPlane * tanf(halfFovy);

    // --- ISOTROPIC DIDACTIC VIEW (default) ---
    float halfWidthNear  = halfHeightNear;
    float halfDepthNdc   = halfHeightNear;
    // // --- ANISOTROPIC VIEW  ---
    // float aspect = (float)screenWidth / (float)screenHeight;
    // float halfWidthNear = halfHeightNear * aspect;
    // float halfDepthNdc  = 0.5f * (farClipPlane - nearClipPlane);

    Vector3 centerNear = Vector3Add(observer->position, Vector3Scale(los, nearClipPlane));
    Vector3 ndcCubeCenter  = Vector3Add(centerNear, Vector3Scale(los, halfDepthNdc));

    drawNdcCube(observer, ndcCubeCenter, halfWidthNear, halfHeightNear, halfDepthNdc, SKYBLUE);

    updateWorldToNdcMappedMesh(&ndcModel->meshes[0], &worldModel->meshes[0],
                               screenWidth, screenHeight, observer, nearClipPlane, farClipPlane,
                               modelPosition, modelScale, meshRotationRadians);
}

//----------------------------------------------------------------------------------
// Topology
//----------------------------------------------------------------------------------
typedef struct WeldedVertex { int id; } WeldedVertex;

typedef struct Topology {
    int triangleCount;
    unsigned short *triangles;                 // indices
    Vector3 *verticesPerTriangle;              // positions per triangle stored as 3 consecutive Vector3s
    int *weldedVertexId;                       // per-vertex welded id
    WeldedVertex *weldedVerticesPerTriangle;   // per-triangle welded vertices
    int (*neighborsPerTriangle)[3];            // neighbor triangle indices per edge
    unsigned char *frontTrianglesFlag;         // bool-like flags for front-facing triangles
    unsigned char *silhouetteTrianglesFlag;    // bool-like flags for silhouette triangles
} Topology;

static void freeTopology(Topology *topology)
{
    if (!topology) return;
    free(topology->triangles);
    free(topology->verticesPerTriangle);
    free(topology->weldedVertexId);
    free(topology->weldedVerticesPerTriangle);
    free(topology->neighborsPerTriangle);
    free(topology->frontTrianglesFlag);
    free(topology->silhouetteTrianglesFlag);
    memset(topology, 0, sizeof(*topology));
}

static void topologyBuild(Topology *topology, Mesh *mesh)
{
    memset(topology, 0, sizeof(topology[0]));
    topology->triangleCount = mesh->triangleCount;
    if (topology->triangleCount <= 0) return;

    topology->triangles = (unsigned short*)malloc(sizeof(unsigned short) * 3 * topology->triangleCount);
    if (mesh->indices) memcpy(topology->triangles, mesh->indices, sizeof(unsigned short) * 3 * topology->triangleCount);
    else for (int linearIndex = 0; linearIndex < 3 * topology->triangleCount; linearIndex++) topology->triangles[linearIndex] = (unsigned short)linearIndex;

    topology->verticesPerTriangle = (Vector3*)malloc(sizeof(Vector3) * 3 * topology->triangleCount);
    for (int triangleIndex = 0; triangleIndex < topology->triangleCount; triangleIndex++)
    {
        unsigned short indexA = topology->triangles[3*triangleIndex+0];
        unsigned short indexB = topology->triangles[3*triangleIndex+1];
        unsigned short indexC = topology->triangles[3*triangleIndex+2];
        topology->verticesPerTriangle[3*triangleIndex+0] = (Vector3){ mesh->vertices[3*indexA+0], mesh->vertices[3*indexA+1], mesh->vertices[3*indexA+2] };
        topology->verticesPerTriangle[3*triangleIndex+1] = (Vector3){ mesh->vertices[3*indexB+0], mesh->vertices[3*indexB+1], mesh->vertices[3*indexB+2] };
        topology->verticesPerTriangle[3*triangleIndex+2] = (Vector3){ mesh->vertices[3*indexC+0], mesh->vertices[3*indexC+1], mesh->vertices[3*indexC+2] };
    }

    topology->weldedVertexId = (int*)malloc(sizeof(int) * mesh->vertexCount);
    for (int vertexIndex = 0; vertexIndex < mesh->vertexCount; vertexIndex++)
    {
        int quantX = (int)lroundf(mesh->vertices[3*vertexIndex+0] / 1e-5f);
        int quantY = (int)lroundf(mesh->vertices[3*vertexIndex+1] / 1e-5f);
        int quantZ = (int)lroundf(mesh->vertices[3*vertexIndex+2] / 1e-5f);
        int foundWelded = -1;
        for (int previous = 0; previous < vertexIndex; previous++)
        {
            int quantXPrev = (int)lroundf(mesh->vertices[3*previous+0] / 1e-5f);
            int quantYPrev = (int)lroundf(mesh->vertices[3*previous+1] / 1e-5f);
            int quantZPrev = (int)lroundf(mesh->vertices[3*previous+2] / 1e-5f);
            if (quantX == quantXPrev && quantY == quantYPrev && quantZ == quantZPrev) { foundWelded = topology->weldedVertexId[previous]; break; }
        }
        topology->weldedVertexId[vertexIndex] = (foundWelded >= 0)? foundWelded : vertexIndex;
    }

    topology->weldedVerticesPerTriangle = (WeldedVertex*)malloc(sizeof(WeldedVertex) * 3 * topology->triangleCount);
    for (int triangleIndex = 0; triangleIndex < topology->triangleCount; triangleIndex++)
    {
        unsigned short indexA = topology->triangles[3*triangleIndex+0];
        unsigned short indexB = topology->triangles[3*triangleIndex+1];
        unsigned short indexC = topology->triangles[3*triangleIndex+2];
        topology->weldedVerticesPerTriangle[3*triangleIndex+0].id = topology->weldedVertexId[indexA];
        topology->weldedVerticesPerTriangle[3*triangleIndex+1].id = topology->weldedVertexId[indexB];
        topology->weldedVerticesPerTriangle[3*triangleIndex+2].id = topology->weldedVertexId[indexC];
    }

    topology->neighborsPerTriangle = (int(*)[3])malloc(sizeof(int) * 3 * topology->triangleCount);
    for (int triangleIndex = 0; triangleIndex < topology->triangleCount; triangleIndex++)
        topology->neighborsPerTriangle[triangleIndex][0] = topology->neighborsPerTriangle[triangleIndex][1] = topology->neighborsPerTriangle[triangleIndex][2] = -1;

    for (int triangleIndex = 0; triangleIndex < topology->triangleCount; triangleIndex++)
    {
        int edge0[2] = { topology->weldedVerticesPerTriangle[3*triangleIndex+0].id, topology->weldedVerticesPerTriangle[3*triangleIndex+1].id };
        int edge1[2] = { topology->weldedVerticesPerTriangle[3*triangleIndex+1].id, topology->weldedVerticesPerTriangle[3*triangleIndex+2].id };
        int edge2[2] = { topology->weldedVerticesPerTriangle[3*triangleIndex+2].id, topology->weldedVerticesPerTriangle[3*triangleIndex+0].id };

        for (int otherTriangleIndex = triangleIndex+1; otherTriangleIndex < topology->triangleCount; otherTriangleIndex++)
        {
            int otherEdge0[2] = { topology->weldedVerticesPerTriangle[3*otherTriangleIndex+0].id, topology->weldedVerticesPerTriangle[3*otherTriangleIndex+1].id };
            int otherEdge1[2] = { topology->weldedVerticesPerTriangle[3*otherTriangleIndex+1].id, topology->weldedVerticesPerTriangle[3*otherTriangleIndex+2].id };
            int otherEdge2[2] = { topology->weldedVerticesPerTriangle[3*otherTriangleIndex+2].id, topology->weldedVerticesPerTriangle[3*otherTriangleIndex+0].id };

            int weldedEdgePairsThis[3][2]  = { {edge0[0],edge0[1]}, {edge1[0],edge1[1]}, {edge2[0],edge2[1]} };
            int weldedEdgePairsOther[3][2] = { {otherEdge0[0],otherEdge0[1]}, {otherEdge1[0],otherEdge1[1]}, {otherEdge2[0],otherEdge2[1]} };

            for (int edgeIndexThis = 0; edgeIndexThis < 3; edgeIndexThis++)
            for (int edgeIndexOther = 0; edgeIndexOther < 3; edgeIndexOther++)
            {
                int edge0VMin = weldedEdgePairsThis[edgeIndexThis][0];
                int edge0VMax = weldedEdgePairsThis[edgeIndexThis][1];
                if (edge0VMin > edge0VMax) { int swapVertex = edge0VMin; edge0VMin = edge0VMax; edge0VMax = swapVertex; }

                int edge1VMin = weldedEdgePairsOther[edgeIndexOther][0];
                int edge1VMax = weldedEdgePairsOther[edgeIndexOther][1];
                if (edge1VMin > edge1VMax) { int swapVertex = edge1VMin; edge1VMin = edge1VMax; edge1VMax = swapVertex; }

                if (edge0VMin == edge1VMin && edge0VMax == edge1VMax) {
                    topology->neighborsPerTriangle[triangleIndex][edgeIndexThis] = otherTriangleIndex;
                    topology->neighborsPerTriangle[otherTriangleIndex][edgeIndexOther] = triangleIndex;
                }
            }
        }
    }

    topology->frontTrianglesFlag      = (unsigned char*)calloc(topology->triangleCount, 1);
    topology->silhouetteTrianglesFlag = (unsigned char*)calloc(topology->triangleCount, 1);
}

static void topologyFrontTriangles(Topology *topology, float rotationRadians, Camera3D *observer)
{
    Vector3 los = observedLineOfSight(observer);
    for (int triangleIndex = 0; triangleIndex < topology->triangleCount; triangleIndex++)
    {
        Vector3 triVerts[3] = {
            topology->verticesPerTriangle[3*triangleIndex+0],
            topology->verticesPerTriangle[3*triangleIndex+1],
            topology->verticesPerTriangle[3*triangleIndex+2]
        };
        rotateVerticesInPlaneSlice(triVerts, 3, rotationRadians);
        Vector3 normal = triangleNormal(triVerts[0], triVerts[1], triVerts[2]);
        if (Vector3DotProduct(normal, los) <= 0.0f) topology->frontTrianglesFlag[triangleIndex] = 1;
    }
}

static void topologySilhouetteTriangles(Topology *topology)
{
    for (int triangleIndex = 0; triangleIndex < topology->triangleCount; triangleIndex++)
    {
        for (int edgeIndex = 0; edgeIndex < 3; edgeIndex++)
        {
            int neighborTriangle = topology->neighborsPerTriangle[triangleIndex][edgeIndex];
            if (neighborTriangle < 0) continue;
            if (topology->frontTrianglesFlag[triangleIndex] != topology->frontTrianglesFlag[neighborTriangle]) {
                topology->silhouetteTrianglesFlag[triangleIndex] = 1;
                topology->silhouetteTrianglesFlag[neighborTriangle] = 1;
            }
        }
    }
}

//----------------------------------------------------------------------------------
// Near-plane intersection “disk” points
//----------------------------------------------------------------------------------
static Vector3 flipYInNearClipPlane(Vector3 intersectionCoord, Camera3D *observer, float nearClipPlane)
{
    Vector3 los, right, up; observedBasis(observer, &los, &right, &up);
    Vector3 centerNear = Vector3Add(observer->position, Vector3Scale(los, nearClipPlane));
    Vector3 toClipPlane = Vector3Subtract(intersectionCoord, centerNear);
    float xComponent = Vector3DotProduct(toClipPlane, right);
    float yComponent = Vector3DotProduct(toClipPlane, up);
    return Vector3Add(centerNear, Vector3Add(Vector3Scale(right, xComponent), Vector3Scale(up, -yComponent)));
}

static void drawNearPlaneIntersectionalDiskMesh(Camera3D *observer,
                                                float nearClipPlane,
                                                Model *nearPlaneIntersectionalDiskModel,
                                                Vector3 modelPosition,
                                                Vector3 modelScale,
                                                float meshRotationRadians,
                                                Topology *topology,
                                                bool reflectYAxis)
{
    Mesh *diskMesh = &nearPlaneIntersectionalDiskModel->meshes[0];
    int diskMeshCapacityVertexCount = diskMesh->vertexCount;

    int intersectionWriteIndex = 0;
    for (int triangleIndex = 0; triangleIndex < topology->triangleCount; triangleIndex++)
    {
        if (!topology->frontTrianglesFlag[triangleIndex]) continue;

        Vector3 triVerts[3] = { topology->verticesPerTriangle[3*triangleIndex+0],
                                topology->verticesPerTriangle[3*triangleIndex+1],
                                topology->verticesPerTriangle[3*triangleIndex+2] };
        rotateVerticesInPlaneSlice(triVerts, 3, meshRotationRadians);

        for (int localVertexIndex = 0; localVertexIndex < 3; localVertexIndex++)
        {
            Vector3 scaledVertex = { triVerts[localVertexIndex].x * modelScale.x, triVerts[localVertexIndex].y * modelScale.y, triVerts[localVertexIndex].z * modelScale.z };
            Vector3 worldVertex = Vector3Add(modelPosition, scaledVertex);
            Vector3 intersectionCoord = nearPlaneIntersection(observer, nearClipPlane, worldVertex);
            if (reflectYAxis) intersectionCoord = flipYInNearClipPlane(intersectionCoord, observer, nearClipPlane);

            rlSetLineWidth(1.0f);
            DrawLine3D(worldVertex, intersectionCoord, (Color){ 255, 0, 0, 80 });

            if (intersectionWriteIndex < diskMeshCapacityVertexCount)
            {
                diskMesh->vertices[3*intersectionWriteIndex+0] = intersectionCoord.x;
                diskMesh->vertices[3*intersectionWriteIndex+1] = intersectionCoord.y;
                diskMesh->vertices[3*intersectionWriteIndex+2] = intersectionCoord.z;
                intersectionWriteIndex++;
            }
        }
    }

    if (intersectionWriteIndex > 0) {
        diskMesh->vertexCount = intersectionWriteIndex;
        UploadMesh(diskMesh, false);
    }
    rlSetPointSize(6.0f);
    DrawModelPoints(*nearPlaneIntersectionalDiskModel, (Vector3){0,0,0}, 1.0f, GREEN);
}

//----------------------------------------------------------------------------------
// Intentionally incorrect projection (affine UV on near plane)
//----------------------------------------------------------------------------------
static void perspectiveIncorrectProjectionDidactic(int screenWidth, int screenHeight,
                                                   Camera3D *observer, float nearClipPlane,
                                                   Mesh *mesh,
                                                   Vector3 modelPosition, Vector3 modelScale,
                                                   float meshRotationRadians, Texture2D texture)
{
    (void)screenWidth; (void)screenHeight; // TODO: not used yet...

    Vector3 los, right, up; observedBasis(observer, &los, &right, &up);

    Vector3 centerNear = Vector3Add(observer->position, Vector3Scale(los, nearClipPlane));
    float halfFovy = DEG2RAD * observer->fovy * 0.5f;
    float halfHeightNear = nearClipPlane * tanf(halfFovy);
    float halfWidthNear  = halfHeightNear;

    rlColor4f(1,1,1,1);
    rlEnableTexture(texture.id);
    rlSetTexture(texture.id);
    rlBegin(RL_TRIANGLES);

    for (int triangleIndex = 0; triangleIndex < mesh->triangleCount; triangleIndex++)
    {
        int indexA = mesh->indices ? mesh->indices[3*triangleIndex+0] : 3*triangleIndex+0;
        int indexB = mesh->indices ? mesh->indices[3*triangleIndex+1] : 3*triangleIndex+1;
        int indexC = mesh->indices ? mesh->indices[3*triangleIndex+2] : 3*triangleIndex+2;

        Vector3 vertexA = { mesh->vertices[3*indexA+0], mesh->vertices[3*indexA+1], mesh->vertices[3*indexA+2] };
        Vector3 vertexB = { mesh->vertices[3*indexB+0], mesh->vertices[3*indexB+1], mesh->vertices[3*indexB+2] };
        Vector3 vertexC = { mesh->vertices[3*indexC+0], mesh->vertices[3*indexC+1], mesh->vertices[3*indexC+2] };

        vertexA = applyModelTranslateRotateScale(vertexA, modelPosition, modelScale, meshRotationRadians);
        vertexB = applyModelTranslateRotateScale(vertexB, modelPosition, modelScale, meshRotationRadians);
        vertexC = applyModelTranslateRotateScale(vertexC, modelPosition, modelScale, meshRotationRadians);

        Vector3 hitA = nearPlaneIntersection(observer, nearClipPlane, vertexA);
        Vector3 hitB = nearPlaneIntersection(observer, nearClipPlane, vertexB);
        Vector3 hitC = nearPlaneIntersection(observer, nearClipPlane, vertexC);

        Vector3 fromNearA = Vector3Subtract(hitA, centerNear);
        Vector3 fromNearB = Vector3Subtract(hitB, centerNear);
        Vector3 fromNearC = Vector3Subtract(hitC, centerNear);

        float sA = 0.5f + 0.5f * Vector3DotProduct(fromNearA, right) / halfWidthNear;
        float tA = 0.5f + 0.5f * Vector3DotProduct(fromNearA, up)    / halfHeightNear;
        float sB = 0.5f + 0.5f * Vector3DotProduct(fromNearB, right) / halfWidthNear;
        float tB = 0.5f + 0.5f * Vector3DotProduct(fromNearB, up)    / halfHeightNear;
        float sC = 0.5f + 0.5f * Vector3DotProduct(fromNearC, right) / halfWidthNear;
        float tC = 0.5f + 0.5f * Vector3DotProduct(fromNearC, up)    / halfHeightNear;

        float uA = mesh->texcoords ? mesh->texcoords[2*indexA+0] : sA;
        float vA = mesh->texcoords ? mesh->texcoords[2*indexA+1] : tA;
        float uB = mesh->texcoords ? mesh->texcoords[2*indexB+0] : sB;
        float vB = mesh->texcoords ? mesh->texcoords[2*indexB+1] : tB;
        float uC = mesh->texcoords ? mesh->texcoords[2*indexC+0] : sC;
        float vC = mesh->texcoords ? mesh->texcoords[2*indexC+1] : tC;

        const float shiftEps = 0.0001f;
        Vector3 planeA = Vector3Add(hitA, Vector3Scale(los, shiftEps));
        Vector3 planeB = Vector3Add(hitB, Vector3Scale(los, shiftEps));
        Vector3 planeC = Vector3Add(hitC, Vector3Scale(los, shiftEps));

        rlTexCoord2f(uA, vA); rlVertex3f(planeA.x, planeA.y, planeA.z);
        rlTexCoord2f(uB, vB); rlVertex3f(planeB.x, planeB.y, planeB.z);
        rlTexCoord2f(uC, vC); rlVertex3f(planeC.x, planeC.y, planeC.z);
    }

    rlEnd();
    rlDisableTexture();
}

//----------------------------------------------------------------------------------
// Palette coloring by triangle; initializes colors and texcoords if needed
//----------------------------------------------------------------------------------
static void applyBarycentricPalette(Mesh *mesh)
{
    if (mesh->colors == NULL)
    {
        mesh->colors = (unsigned char*)RL_CALLOC(mesh->vertexCount * 4, sizeof(unsigned char));
    }
    for (int triangleIndex = 0; triangleIndex < mesh->triangleCount; triangleIndex++)
    {
        int indexA = mesh->indices ? mesh->indices[3*triangleIndex+0] : 3*triangleIndex+0;
        int indexB = mesh->indices ? mesh->indices[3*triangleIndex+1] : 3*triangleIndex+1;
        int indexC = mesh->indices ? mesh->indices[3*triangleIndex+2] : 3*triangleIndex+2;

        mesh->colors[4*indexA+0] = 255; mesh->colors[4*indexA+1] = 0;   mesh->colors[4*indexA+2] = 0;   mesh->colors[4*indexA+3] = 255;
        mesh->colors[4*indexB+0] = 0;   mesh->colors[4*indexB+1] = 255; mesh->colors[4*indexB+2] = 0;   mesh->colors[4*indexB+3] = 255;
        mesh->colors[4*indexC+0] = 0;   mesh->colors[4*indexC+1] = 0;   mesh->colors[4*indexC+2] = 255; mesh->colors[4*indexC+3] = 255;
    }

    if (mesh->texcoords == NULL)
    {
        mesh->texcoords = (float*)RL_CALLOC(mesh->vertexCount * 2, sizeof(float));
    }
    for (int triangleIndex = 0; triangleIndex < mesh->triangleCount; triangleIndex++)
    {
        int indexA = mesh->indices ? mesh->indices[3*triangleIndex+0] : 3*triangleIndex+0;
        int indexB = mesh->indices ? mesh->indices[3*triangleIndex+1] : 3*triangleIndex+1;
        int indexC = mesh->indices ? mesh->indices[3*triangleIndex+2] : 3*triangleIndex+2;

        mesh->texcoords[2*indexA+0] = 1.0f; mesh->texcoords[2*indexA+1] = 0.0f;
        mesh->texcoords[2*indexB+0] = 0.0f; mesh->texcoords[2*indexB+1] = 1.0f;
        mesh->texcoords[2*indexC+0] = 0.0f; mesh->texcoords[2*indexC+1] = 0.0f;
    }

    UploadMesh(mesh, false);
}

static void drawQuad(Vector3 a, Vector3 b, Vector3 c, Vector3 d, Color color)
{
    DrawTriangle3D(a, b, c, color);
    DrawTriangle3D(a, c, d, color);
}

// ---------- Software rasterizer visualization (didactic, draws near-plane “pixels” in 3D) ----------
static Vector3 perspectiveProjectStageData(int screenWidth, int screenHeight,
                                           Camera3D *observer, float nearClipPlane,
                                           Vector3 worldCoord)
{
    float aspect = (float)screenWidth/(float)screenHeight;
    Vector3 los, right, up; observedBasis(observer, &los, &right, &up);

    Vector3 centerNear = Vector3Add(observer->position, Vector3Scale(los, nearClipPlane));
    float halfFovy = DEG2RAD * observer->fovy * 0.5f;
    float halfHeightNear = nearClipPlane * tanf(halfFovy);
    float halfWidthNear  = halfHeightNear * aspect;

    Vector3 rayToWorld = Vector3Subtract(worldCoord, observer->position);
    float signedDepth = Vector3DotProduct(rayToWorld, los);

    float xNdc, yNdc;
#if PERSPECTIVE_CORRECT
    {
        Vector3 intersectionCoord = nearPlaneIntersection(observer, nearClipPlane, worldCoord);
        Vector3 clipPlaneVector = Vector3Subtract(intersectionCoord, centerNear);
        xNdc = Vector3DotProduct(clipPlaneVector, right) / (halfWidthNear + 1e-9f);
        yNdc = Vector3DotProduct(clipPlaneVector, up)    / (halfHeightNear + 1e-9f);
    }
#else
    {
        Vector3 dir = Vector3Normalize(rayToWorld);
        xNdc = Vector3DotProduct(dir, right) / (halfWidthNear + 1e-9f);
        yNdc = Vector3DotProduct(dir, up)    / (halfHeightNear + 1e-9f);
    }
#endif

    float pixelX = (xNdc*0.5f + 0.5f)*(float)screenWidth;
    float pixelY = (1.0f - (yNdc*0.5f + 0.5f))*(float)screenHeight;
    return (Vector3){ pixelX, pixelY, signedDepth };
}

static inline float fullSignedBarycentricWeightComponent(Vector2 A, Vector2 B, Vector2 C)
{
    return (C.x - A.x)*(B.y - A.y) - (C.y - A.y)*(B.x - A.x);
}

static inline float partialSignedBarycentricWeightComponent(Vector2 edgeStart, Vector2 edgeEnd, Vector2 sample)
{
    return (sample.x - edgeStart.x)*(edgeEnd.y - edgeStart.y) - (sample.y - edgeStart.y)*(edgeEnd.x - edgeStart.x);
}

static inline bool barycentricWeightsForSample(Vector2 sample, Vector2 A, Vector2 B, Vector2 C,
                                               float denom, float *weightsOut)
{
    float wA = partialSignedBarycentricWeightComponent(B, C, sample);
    float wB = partialSignedBarycentricWeightComponent(C, A, sample);
    float wC = partialSignedBarycentricWeightComponent(A, B, sample);

    if ((wA*denom >= 0.0f) && (wB*denom >= 0.0f) && (wC*denom >= 0.0f))
    {
        weightsOut[0] = wA/denom; weightsOut[1] = wB/denom; weightsOut[2] = wC/denom;
        return true;
    }
    return false;
}

static inline float interpolateDepth(const float *weights, const float *inverseDepths)
{
    float sum = weights[0]*inverseDepths[0] + weights[1]*inverseDepths[1] + weights[2]*inverseDepths[2];
    return 1.0f/(sum + 1e-9f);
}

static inline Color interpolateColor3(const float *weights, int indexA, int indexB, int indexC, const unsigned char *colors)
{
    float r = weights[0]*colors[4*indexA+0] + weights[1]*colors[4*indexB+0] + weights[2]*colors[4*indexC+0];
    float g = weights[0]*colors[4*indexA+1] + weights[1]*colors[4*indexB+1] + weights[2]*colors[4*indexC+1];
    float b = weights[0]*colors[4*indexA+2] + weights[1]*colors[4*indexB+2] + weights[2]*colors[4*indexC+2];
    return (Color){ (unsigned char)Clamp(r,0,255), (unsigned char)Clamp(g,0,255), (unsigned char)Clamp(b,0,255), 255 };
}

static inline void rasterizationBoundingBox(Vector2 A, Vector2 B, Vector2 C, int screenWidth, int screenHeight,
                                            int *minXOut, int *maxXOut, int *minYOut, int *maxYOut)
{
    float minx = fminf(A.x, fminf(B.x, C.x));
    float maxx = fmaxf(A.x, fmaxf(B.x, C.x));
    float miny = fminf(A.y, fminf(B.y, C.y));
    float maxy = fmaxf(A.y, fmaxf(B.y, C.y));

    int ix0 = (int)floorf(minx), ix1 = (int)ceilf(maxx);
    int iy0 = (int)floorf(miny), iy1 = (int)ceilf(maxy);

    if (ix0 < 0) ix0 = 0; if (iy0 < 0) iy0 = 0;
    if (ix1 > screenWidth-1) ix1 = screenWidth-1;
    if (iy1 > screenHeight-1) iy1 = screenHeight-1;

    *minXOut = ix0; *maxXOut = ix1; *minYOut = iy0; *maxYOut = iy1;
}

static inline bool depthTestWrite(float depth, float *depthBuffer, int bufferIndex)
{
#if DEPTH_TEST_ON
    if (depth < depthBuffer[bufferIndex])
    {
        depthBuffer[bufferIndex] = depth;
        return true;
    }
    return false;
#else
    (void)depthBuffer; (void)bufferIndex; (void)depth;
    return true;
#endif
}

static void rasterizeTriangle(Camera3D *observer, int screenWidth, int screenHeight, float nearClipPlane,
                              Vector3 worldA, Vector3 worldB, Vector3 worldC,
                              int indexA, int indexB, int indexC,
                              Mesh *mesh, float *depthBuffer, int step)
{
    Vector3 projA = perspectiveProjectStageData(screenWidth, screenHeight, observer, nearClipPlane, worldA);
    Vector3 projB = perspectiveProjectStageData(screenWidth, screenHeight, observer, nearClipPlane, worldB);
    Vector3 projC = perspectiveProjectStageData(screenWidth, screenHeight, observer, nearClipPlane, worldC);

    Vector2 A = { projA.x, projA.y }, B = { projB.x, projB.y }, C = { projC.x, projC.y };
    float z[3] = { projA.z, projB.z, projC.z };
    float invZ[3] = { 1.0f/(z[0]+1e-9f), 1.0f/(z[1]+1e-9f), 1.0f/(z[2]+1e-9f) };

    float baryDenom = fullSignedBarycentricWeightComponent(A, B, C);
    if (fabsf(baryDenom) <= 1e-9f) return;

    int minX, maxX, minY, maxY;
    rasterizationBoundingBox(A, B, C, screenWidth, screenHeight, &minX, &maxX, &minY, &maxY);

    Vector3 los, right, up; observedBasis(observer, &los, &right, &up);

    float halfFovy = DEG2RAD * observer->fovy * 0.5f;
    float halfHeightNear = nearClipPlane * tanf(halfFovy);
    float halfWidthNear  = halfHeightNear * ((float)screenWidth / (float)screenHeight);

    Vector3 centerNear = Vector3Add(observer->position, Vector3Scale(los, nearClipPlane));
    Vector3 pixelStepHalfX = Vector3Scale(right, (halfWidthNear / (float)screenWidth));
    Vector3 pixelStepHalfY = Vector3Scale(up,    (halfHeightNear / (float)screenHeight));

    for (int pixelY = (minY/step)*step; pixelY <= maxY; pixelY += step)
    {
        float sampleY = (float)pixelY + 0.5f;
        for (int pixelX = (minX/step)*step; pixelX <= maxX; pixelX += step)
        {
            float sampleX = (float)pixelX + 0.5f;
            float weights[3];
            if (barycentricWeightsForSample((Vector2){sampleX, sampleY}, A, B, C, baryDenom, weights))
            {
                float depth = interpolateDepth(weights, invZ);
                int bufferIndex = pixelY * screenWidth + pixelX;
                if (depthTestWrite(depth, depthBuffer, bufferIndex))
                {
                    Color color = (Color){ 255,255,255,255 };
                    if (mesh->colors)
                    {
                        color = interpolateColor3(weights, indexA, indexB, indexC, mesh->colors);
                        color.a = 255;
                    }

                    float xNdc = (sampleX/(float)screenWidth)*2.0f - 1.0f;
                    float yNdc = 1.0f - (sampleY/(float)screenHeight)*2.0f;

                    Vector3 nearOrigin = Vector3Add(centerNear,
                                         Vector3Add(Vector3Scale(right, xNdc * halfWidthNear),
                                                    Vector3Scale(up,    yNdc * halfHeightNear)));

                    Vector3 quadCornerA = Vector3Subtract(Vector3Subtract(nearOrigin, pixelStepHalfX), pixelStepHalfY);
                    Vector3 quadCornerB = Vector3Add(     Vector3Subtract(nearOrigin, pixelStepHalfY), pixelStepHalfX);
                    Vector3 quadCornerC = Vector3Add(     Vector3Add(nearOrigin, pixelStepHalfX),     pixelStepHalfY);
                    Vector3 quadCornerD = Vector3Add(     Vector3Subtract(nearOrigin, pixelStepHalfX), pixelStepHalfY);

                    DrawTriangle3D(quadCornerA, quadCornerB, quadCornerC, color);
                    DrawTriangle3D(quadCornerA, quadCornerC, quadCornerD, color);
                }
            }
        }
    }
}

static void drawNearPlaneSoftwareRaster(Camera3D *observer,
                                        int screenWidth, int screenHeight,
                                        float nearClipPlane,
                                        Mesh *mesh,
                                        Vector3 modelPosition, Vector3 modelScale,
                                        float meshRotationRadians,
                                        int step)
{
    float *depthBuffer = (float*)RL_CALLOC((size_t)screenWidth * screenHeight, sizeof(float));
    for (int i = 0; i < screenWidth * screenHeight; i++) depthBuffer[i] = INFINITY;

    float s = sinf(meshRotationRadians), c = cosf(meshRotationRadians);

    for (int triangleIndex = 0; triangleIndex < mesh->triangleCount; triangleIndex++)
    {
        int indexA = mesh->indices ? mesh->indices[3*triangleIndex+0] : 3*triangleIndex+0;
        int indexB = mesh->indices ? mesh->indices[3*triangleIndex+1] : 3*triangleIndex+1;
        int indexC = mesh->indices ? mesh->indices[3*triangleIndex+2] : 3*triangleIndex+2;

        Vector3 vertexA = { mesh->vertices[3*indexA+0], mesh->vertices[3*indexA+1], mesh->vertices[3*indexA+2] };
        Vector3 vertexB = { mesh->vertices[3*indexB+0], mesh->vertices[3*indexB+1], mesh->vertices[3*indexB+2] };
        Vector3 vertexC = { mesh->vertices[3*indexC+0], mesh->vertices[3*indexC+1], mesh->vertices[3*indexC+2] };

        Vector3 tri[3] = { vertexA, vertexB, vertexC };
        for (int k = 0; k < 3; k++)
        {
            float x0 = tri[k].x, z0 = tri[k].z;
            tri[k].x = c * x0 + s * z0;
            tri[k].z = -s * x0 + c * z0;
        }

        for (int k = 0; k < 3; k++)
        {
            tri[k].x = tri[k].x * modelScale.x + modelPosition.x;
            tri[k].y = tri[k].y * modelScale.y + modelPosition.y;
            tri[k].z = tri[k].z * modelScale.z + modelPosition.z;
        }

        rasterizeTriangle(observer, screenWidth, screenHeight, nearClipPlane,
                          tri[0], tri[1], tri[2],
                          indexA, indexB, indexC, mesh, depthBuffer, step);
    }

    RL_FREE(depthBuffer);
}

//----------------------------------------------------------------------------------
// Jugemu orbital camera
//----------------------------------------------------------------------------------
static bool moveJugemuOrbital(Camera3D *jugemu, float deltaTime)
{
    Vector3 previousPosition = jugemu->position;

    float radius = Vector3Length(jugemu->position);
    float azimuth = atan2f(jugemu->position.z, jugemu->position.x);
    float horizontalRadius = sqrtf(jugemu->position.x * jugemu->position.x + jugemu->position.z * jugemu->position.z);
    float elevation = atan2f(jugemu->position.y, horizontalRadius);

    const float LONG_SPEED = 1.5f;
    const float LAT_SPEED  = 1.0f;
    const float ZOOM_SPEED = 2.0f;
    const float ROLL_SPEED = 2.0f;

    if (IsKeyDown(KEY_LEFT))  azimuth   += LONG_SPEED*deltaTime;
    if (IsKeyDown(KEY_RIGHT)) azimuth   -= LONG_SPEED*deltaTime;
    if (IsKeyDown(KEY_UP))    elevation += LAT_SPEED*deltaTime;
    if (IsKeyDown(KEY_DOWN))  elevation -= LAT_SPEED*deltaTime;
    if (IsKeyDown(KEY_W))     radius    -= ZOOM_SPEED*deltaTime;
    if (IsKeyDown(KEY_S))     radius    += ZOOM_SPEED*deltaTime;

    float rollDeltaRadians = 0.0f;
    if (IsKeyDown(KEY_A)) rollDeltaRadians -= ROLL_SPEED*deltaTime;
    if (IsKeyDown(KEY_D)) rollDeltaRadians += ROLL_SPEED*deltaTime;
    if (IsKeyPressed(KEY_SPACE)) { jugemu->up = (Vector3){0,1,0}; rollDeltaRadians = 0.0f; }

    radius = Clamp(radius, 0.25f, 25.0f);
    const float EPS = 0.0001f;
    elevation = Clamp(elevation, -PI*0.5f + EPS, PI*0.5f - EPS);

    jugemu->position.x = radius * cosf(elevation) * cosf(azimuth);
    jugemu->position.y = radius * sinf(elevation);
    jugemu->position.z = radius * cosf(elevation) * sinf(azimuth);

    Vector3 viewDir = Vector3Normalize(Vector3Subtract((Vector3){0,0,0}, jugemu->position));
    Vector3 upRot = rotatePointAboutAxis(jugemu->up, (Vector3){0,0,0}, viewDir, rollDeltaRadians);
    jugemu->target = (Vector3){0,0,0};
    jugemu->up = Vector3Normalize(upRot);

    return (previousPosition.x != jugemu->position.x) || (previousPosition.y != jugemu->position.y) || (previousPosition.z != jugemu->position.z);
}

//------------------------------------------------------------------------------------
// Program main entry point
//------------------------------------------------------------------------------------
int main(void)
{
    // Initialization
    //--------------------------------------------------------------------------------------
    InitWindow(N64_WIDTH, N64_HEIGHT, "raylib teaching - world->NDC pipeline (single-file)");
    SetTargetFPS(60);

    float nearClipPlane = 1.0f;
    float farClipPlane  = 3.0f;

    Camera3D mainObserver = { 0 };
    mainObserver.position = OBSERVER_POS;
    mainObserver.target   = (Vector3){ 0, 0, 0 };
    mainObserver.up       = (Vector3){ 0, 1, 0 };
    mainObserver.fovy     = FOVY_PERSPECTIVE;
    mainObserver.projection = CAMERA_PERSPECTIVE;

    int screenWidth  = GetScreenWidth();
    int screenHeight = GetScreenHeight();
    float aspect = (float)screenWidth/(float)screenHeight;

    Camera3D jugemu = (Camera3D){ 0 };
    jugemu.position = JUGEMU_POS_ISO;
    jugemu.target   = (Vector3){ 0, 0, 0 };
    jugemu.up       = (Vector3){ 0, 1, 0 };
    jugemu.fovy     = FOVY_PERSPECTIVE;
    jugemu.projection = CAMERA_PERSPECTIVE;

    TraceLog(LOG_INFO, TextFormat("jugemu init pos: (%.3f, %.3f, %.3f)", jugemu.position.x, jugemu.position.y, jugemu.position.z));

    float idleTimer = 0.0f;
    bool movedSinceLastLog = false;

    float meshRotationRadians = 0.0f;

    Mesh cubeMesh = GenMeshCube(1.0f, 1.0f, 1.0f);
    Model mainModel = LoadModelFromMesh(cubeMesh);

    Mesh ndcMesh = (Mesh){ 0 };
    ndcMesh.vertexCount   = mainModel.meshes[0].vertexCount;
    ndcMesh.triangleCount = mainModel.meshes[0].triangleCount;
    ndcMesh.vertices  = RL_CALLOC(ndcMesh.vertexCount*3, sizeof(float));
    ndcMesh.texcoords = RL_CALLOC(ndcMesh.vertexCount*2, sizeof(float));
    ndcMesh.indices   = RL_CALLOC(ndcMesh.triangleCount*3, sizeof(unsigned short));
    if (mainModel.meshes[0].texcoords)
    {
        for (int vertexIndex = 0; vertexIndex < ndcMesh.vertexCount; vertexIndex++)
        {
            ndcMesh.texcoords[2*vertexIndex+0] = mainModel.meshes[0].texcoords[2*vertexIndex+0];
            ndcMesh.texcoords[2*vertexIndex+1] = mainModel.meshes[0].texcoords[2*vertexIndex+1];
        }
    }
    UploadMesh(&ndcMesh, false);
    Model ndcModel = LoadModelFromMesh(ndcMesh);

    Mesh nearPlaneIntersectionalDiskMesh = (Mesh){ 0 };
    nearPlaneIntersectionalDiskMesh.vertexCount = mainModel.meshes[0].vertexCount;
    nearPlaneIntersectionalDiskMesh.vertices = RL_CALLOC(nearPlaneIntersectionalDiskMesh.vertexCount*3, sizeof(float));
    UploadMesh(&nearPlaneIntersectionalDiskMesh, false);
    Model nearPlaneIntersectionalDiskModel = LoadModelFromMesh(nearPlaneIntersectionalDiskMesh);

    Image checked = GenImageChecked(16, 16, 4, 4, BLACK, WHITE);
    Texture2D tex = LoadTextureFromImage(checked);
    UnloadImage(checked);
    mainModel.materials[0].maps[MATERIAL_MAP_ALBEDO].texture = tex;
    ndcModel.materials[0].maps[MATERIAL_MAP_ALBEDO].texture  = tex;

    unsigned int NDC_MODE = 0;
    static unsigned int ASPECT_MODE = 0; // start isotropic for teaching
    //--------------------------------------------------------------------------------------

    // Main loop
    while (!WindowShouldClose())
    {
        // Update
        //----------------------------------------------------------------------------------
        if (IsKeyPressed(KEY_E)) NDC_MODE = NDC_MODE + 1;
        if (IsKeyPressed(KEY_Q)) ASPECT_MODE ^= 1u;
        screenWidth  = GetScreenWidth();
        screenHeight = GetScreenHeight();
        aspect       = (float)screenWidth/(float)screenHeight;

        float deltaTime = GetFrameTime();
        meshRotationRadians -= ANGULAR_VELOCITY * deltaTime;

        if (moveJugemuOrbital(&jugemu, deltaTime))
        {
            movedSinceLastLog = true;
            idleTimer = 0.0f;
        }
        else if (movedSinceLastLog)
        {
            idleTimer += deltaTime;
            if (idleTimer >= 1.0f)
            {
                TraceLog(LOG_INFO, TextFormat("jugemu stopped at: (%.3f, %.3f, %.3f)",
                        jugemu.position.x, jugemu.position.y, jugemu.position.z));
                movedSinceLastLog = false;
                idleTimer = 0.0f;
            }
        }
        //----------------------------------------------------------------------------------

        // Draw
        //----------------------------------------------------------------------------------
        BeginDrawing();
            ClearBackground(BLACK);

            BeginMode3D(jugemu);

                drawObservedAxes(&mainObserver);

                if ((NDC_MODE & 1u) == 1u)
                {
                    rlDisableBackfaceCulling();

                    mapFrustumToNdcCube(&mainObserver, screenWidth, screenHeight,
                                        nearClipPlane, farClipPlane,
                                        &mainModel, &ndcModel,
                                        MODEL_POS, MODEL_SCALE, meshRotationRadians);
                    ndcModel.materials[0].maps[MATERIAL_MAP_ALBEDO].texture.id  = tex.id;
                    DrawModelEx(ndcModel, MODEL_POS, (Vector3){0,1,0}, RAD2DEG*meshRotationRadians, MODEL_SCALE, WHITE);
                    rlSetLineWidth(2.0f);
                    ndcModel.materials[0].maps[MATERIAL_MAP_ALBEDO].texture.id  = 0;
                    DrawModelWiresEx(ndcModel, MODEL_POS, (Vector3){0,1,0}, RAD2DEG*meshRotationRadians, MODEL_SCALE, BLUE);
                    rlSetPointSize(8.0f);
                    DrawModelPointsEx(ndcModel, MODEL_POS, (Vector3){0,1,0}, RAD2DEG*meshRotationRadians, MODEL_SCALE, GREEN);

                    Topology topoNdc; topologyBuild(&topoNdc, &ndcModel.meshes[0]);
                    topologyFrontTriangles(&topoNdc, meshRotationRadians, &mainObserver);
                    topologySilhouetteTriangles(&topoNdc);

                    drawNearPlaneIntersectionalDiskMesh(&mainObserver, nearClipPlane,
                        &nearPlaneIntersectionalDiskModel,
                        MODEL_POS, MODEL_SCALE, meshRotationRadians,
                        &topoNdc, true);

                    perspectiveIncorrectProjectionDidactic(screenWidth, screenHeight, &mainObserver,
                        nearClipPlane, &ndcModel.meshes[0], MODEL_POS, MODEL_SCALE, meshRotationRadians, tex);

                    freeTopology(&topoNdc);
                }
                else
                {
                    mainModel.materials[0].maps[MATERIAL_MAP_ALBEDO].texture.id = tex.id;
                    DrawModelEx(mainModel, MODEL_POS, (Vector3){0,1,0}, RAD2DEG*meshRotationRadians, MODEL_SCALE, WHITE);
                    rlSetLineWidth(2.0f);
                    mainModel.materials[0].maps[MATERIAL_MAP_ALBEDO].texture.id = 0;
                    DrawModelWiresEx(mainModel, MODEL_POS, (Vector3){0,1,0}, RAD2DEG*meshRotationRadians, MODEL_SCALE, BLUE);
                    rlSetPointSize(8.0f);
                    DrawModelPointsEx(mainModel, MODEL_POS, (Vector3){0,1,0}, RAD2DEG*meshRotationRadians, MODEL_SCALE, GREEN);

                    drawFrustum(&mainObserver, aspect, nearClipPlane, farClipPlane);

                    Topology topo; topologyBuild(&topo, &mainModel.meshes[0]);
                    topologyFrontTriangles(&topo, meshRotationRadians, &mainObserver);
                    topologySilhouetteTriangles(&topo);

                    drawNearPlaneIntersectionalDiskMesh(&mainObserver, nearClipPlane,
                        &nearPlaneIntersectionalDiskModel,
                        MODEL_POS, MODEL_SCALE, meshRotationRadians,
                        &topo, false);

                    perspectiveIncorrectProjectionDidactic(screenWidth, screenHeight, &mainObserver,
                        nearClipPlane, &mainModel.meshes[0], MODEL_POS, MODEL_SCALE, meshRotationRadians, tex);

                    freeTopology(&topo);
                }

                // Optional software raster
                // int rasterStep = 4;
                // drawNearPlaneSoftwareRaster(&mainObserver, screenWidth, screenHeight, nearClipPlane,
                //                             &mainModel.meshes[0], MODEL_POS, MODEL_SCALE, meshRotationRadians, rasterStep);

            EndMode3D();

            DrawText("E: toggle NDC overlay · Arrows: orbit · W/S: zoom · A/D: roll · Space: reset",
                     12, 14, 10, RAYWHITE);
            DrawText(TextFormat("NDC_MODE bit: %u", (NDC_MODE & 1u)), 12, 28, 10, RAYWHITE);

        EndDrawing();
        //----------------------------------------------------------------------------------
    }

    // De-Initialization
    //--------------------------------------------------------------------------------------
    UnloadModel(mainModel);
    UnloadModel(ndcModel);
    UnloadModel(nearPlaneIntersectionalDiskModel);
    UnloadTexture(tex);
    CloseWindow();
    //--------------------------------------------------------------------------------------

    return 0;
}
