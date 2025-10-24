/*
*******************************************************************************************
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

static const Vector3 MODEL_POS   = {0.0f, 0.0f, 0.0f};
static const Vector3 MODEL_SCALE = {1.0f, 1.0f, 1.0f};

static const Vector3 OBSERVER_POS   = {0.0f, 0.0f, 2.0f};
static const Vector3 JUGEMU_POS_ISO = {3.0f, 1.0f, 3.0f};

//----------------------------------------------------------------------------------
// Global flags/state
//----------------------------------------------------------------------------------
enum
{
    FLAG_NDC_OVERLAY         = 1u << 0, // E
    FLAG_ASPECT              = 1u << 1, // Q  (0=ISO didactic, 1=ANISO true aspect)
    FLAG_PERSPECTIVE_CORRECT = 1u << 2, // P
    FLAG_PAUSE               = 1u << 3, // F
    FLAG_COLOR_MODE          = 1u << 4  // C  <-- NEW: toggle color vs texture
};

static unsigned int gFlags = 0;

#define ANISOTROPIC() ((gFlags & FLAG_ASPECT) != 0)
#define PERSPECTIVE_CORRECT() ((gFlags & FLAG_PERSPECTIVE_CORRECT) != 0)
#define PAUSED() ((gFlags & FLAG_PAUSE) != 0)
#define COLOR_MODE() ((gFlags & FLAG_COLOR_MODE) != 0)
//----------------------------------------------------------------------------------
// "PerspectiveCorrect" screenshot resources CPU capture OpenGL11 + near-plane quad
//----------------------------------------------------------------------------------
static Texture2D gPerspectiveCorrectTexture = {0};
static bool gPerspectiveCorrectTextureReady = false;

static Mesh gNearQuad       = {0};
static Model gNearQuadModel = {0};
static bool gNearQuadBuilt  = false;

//----------------------------------------------------------------------------------
// Topology
//----------------------------------------------------------------------------------
typedef struct WeldedVertex
{
    int id;
} WeldedVertex;

typedef struct Topology
{
    int triangleCount;
    unsigned short *triangles;               // indices
    Vector3 *verticesPerTriangle;            // positions per triangle stored as 3 consecutive Vector3s
    int *weldedVertexId;                     // per-vertex welded id
    WeldedVertex *weldedVerticesPerTriangle; // per-triangle welded vertices
    int (*neighborsPerTriangle)[3];          // neighbor triangle indices per edge
    unsigned char *frontTrianglesFlag;       // bool-like flags for front-facing triangles
    unsigned char *silhouetteTrianglesFlag;  // bool-like flags for silhouette triangles
} Topology;

static bool moveJugemuOrbital(Camera3D *jugemu, float deltaTime);
static void capturePerspectiveCorrectToTexture(
    Camera3D *mainObserver,
    Model *model,
    Vector3 modelPos,
    float rotRadians,
    Vector3 modelScale,
    Texture2D tex,
    bool usePVC);
static void drawObservedAxes(Camera3D *observer);
static void mapFrustumToNdcCube(
    Camera3D *observer,
    int screenWidth,
    int screenHeight,
    float nearClipPlane,
    float farClipPlane,
    Model *worldModel,
    Model *ndcModel,
    Vector3 modelPosition,
    Vector3 modelScale,
    float meshRotationRadians);
static void drawFrustum(Camera3D *observer, float aspect, float nearClipPlane, float farClipPlane);
static void topologyBuild(Topology *topology, Mesh *mesh);
static void topologyFrontTriangles(Topology *topology, float rotationRadians, Camera3D *observer);
static void topologySilhouetteTriangles(Topology *topology);
static void drawNearPlaneIntersectionalDiskMesh(
    Camera3D *observer,
    float nearClipPlane,
    Model *nearPlaneIntersectionalDiskModel,
    Vector3 modelPosition,
    Vector3 modelScale,
    float meshRotationRadians,
    Topology *topology,
    bool reflectYAxis);
static void ensureNearPlaneQuadBuilt(void);
static void updateNearPlaneQuadGeometry(Camera3D *observer, float nearClipPlane, int screenWidth, int screenHeight);
static void applyBarycentricPalette(Mesh *mesh);
static void perspectiveIncorrectColorDidactic(
    int screenWidth,
    int screenHeight,
    Camera3D *observer,
    float nearClipPlane,
    Mesh *mesh,
    Vector3 modelPosition,
    Vector3 modelScale,
    float meshRotationRadians);
static void perspectiveIncorrectProjectionDidactic(
    int screenWidth,
    int screenHeight,
    Camera3D *observer,
    float nearClipPlane,
    Mesh *mesh,
    Vector3 modelPosition,
    Vector3 modelScale,
    float meshRotationRadians,
    Texture2D texture);
static void freeTopology(Topology *topology);

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

    Camera3D mainObserver   = {0};
    mainObserver.position   = OBSERVER_POS;
    mainObserver.target     = (Vector3){0, 0, 0};
    mainObserver.up         = (Vector3){0, 1, 0};
    mainObserver.fovy       = FOVY_PERSPECTIVE;
    mainObserver.projection = CAMERA_PERSPECTIVE;

    int screenWidth  = GetScreenWidth();
    int screenHeight = GetScreenHeight();
    float aspect     = (float)screenWidth / (float)screenHeight;

    Camera3D jugemu   = (Camera3D){0};
    jugemu.position   = JUGEMU_POS_ISO;
    jugemu.target     = (Vector3){0, 0, 0};
    jugemu.up         = (Vector3){0, 1, 0};
    jugemu.fovy       = FOVY_PERSPECTIVE;
    jugemu.projection = CAMERA_PERSPECTIVE;

    TraceLog(
        LOG_INFO,
        TextFormat("jugemu init pos: (%.3f, %.3f, %.3f)", jugemu.position.x, jugemu.position.y, jugemu.position.z));

    float idleTimer        = 0.0f;
    bool movedSinceLastLog = false;

    float meshRotationRadians = 0.0f;

    Mesh cubeMesh   = GenMeshCube(1.0f, 1.0f, 1.0f);
    Model mainModel = LoadModelFromMesh(cubeMesh);

    Mesh ndcMesh          = (Mesh){0};
    ndcMesh.vertexCount   = mainModel.meshes[0].vertexCount;
    ndcMesh.triangleCount = mainModel.meshes[0].triangleCount;
    ndcMesh.vertices      = RL_CALLOC(ndcMesh.vertexCount * 3, sizeof(float));
    ndcMesh.texcoords     = RL_CALLOC(ndcMesh.vertexCount * 2, sizeof(float));
    ndcMesh.indices       = RL_CALLOC(ndcMesh.triangleCount * 3, sizeof(unsigned short));
    if (mainModel.meshes[0].texcoords)
    {
        for (int vertexIndex = 0; vertexIndex < ndcMesh.vertexCount; vertexIndex++)
        {
            ndcMesh.texcoords[2 * vertexIndex + 0] = mainModel.meshes[0].texcoords[2 * vertexIndex + 0];
            ndcMesh.texcoords[2 * vertexIndex + 1] = mainModel.meshes[0].texcoords[2 * vertexIndex + 1];
        }
    }
    UploadMesh(&ndcMesh, false);
    Model ndcModel = LoadModelFromMesh(ndcMesh);

    Mesh nearPlaneIntersectionalDiskMesh        = (Mesh){0};
    nearPlaneIntersectionalDiskMesh.vertexCount = mainModel.meshes[0].vertexCount;
    nearPlaneIntersectionalDiskMesh.vertices =
        RL_CALLOC(nearPlaneIntersectionalDiskMesh.vertexCount * 3, sizeof(float));
    UploadMesh(&nearPlaneIntersectionalDiskMesh, false);
    Model nearPlaneIntersectionalDiskModel = LoadModelFromMesh(nearPlaneIntersectionalDiskMesh);

    Image checked = GenImageChecked(16, 16, 4, 4, BLACK, WHITE);
    Texture2D tex = LoadTextureFromImage(checked);
    UnloadImage(checked);
    mainModel.materials[0].maps[MATERIAL_MAP_ALBEDO].texture = tex;
    ndcModel.materials[0].maps[MATERIAL_MAP_ALBEDO].texture  = tex;
    applyBarycentricPalette(&mainModel.meshes[0]);
    //--------------------------------------------------------------------------------------

    // Main loop
    while (!WindowShouldClose())
    {
        // Update
        //----------------------------------------------------------------------------------
        // Update
        //----------------------------------------------------------------------------------
        if (IsKeyPressed(KEY_E))
            gFlags ^= FLAG_NDC_OVERLAY;
        if (IsKeyPressed(KEY_Q))
            gFlags ^= FLAG_ASPECT;
        if (IsKeyPressed(KEY_P))
            gFlags ^= FLAG_PERSPECTIVE_CORRECT; // <- no auto-freeze, ever
        if (IsKeyPressed(KEY_F))
            gFlags ^= FLAG_PAUSE; // <- only way to pause
        if (IsKeyPressed(KEY_C))
            gFlags ^= FLAG_COLOR_MODE;

        screenWidth  = GetScreenWidth();
        screenHeight = GetScreenHeight();
        aspect       = (float)screenWidth / (float)screenHeight;

        float deltaTime = GetFrameTime();

        if (!PAUSED())
        {
            meshRotationRadians -= ANGULAR_VELOCITY * deltaTime;

            if (moveJugemuOrbital(&jugemu, deltaTime))
            {
                movedSinceLastLog = true;
                idleTimer         = 0.0f;
            }
            else if (movedSinceLastLog)
            {
                idleTimer += deltaTime;
                if (idleTimer >= 1.0f)
                {
                    TraceLog(
                        LOG_INFO,
                        TextFormat(
                            "jugemu stopped at: (%.3f, %.3f, %.3f)",
                            jugemu.position.x,
                            jugemu.position.y,
                            jugemu.position.z));
                    movedSinceLastLog = false;
                    idleTimer         = 0.0f;
                }
            }
        }
        //----------------------------------------------------------------------------------
        if (moveJugemuOrbital(&jugemu, deltaTime))
        {
            movedSinceLastLog = true;
            idleTimer         = 0.0f;
        }
        else if (movedSinceLastLog)
        {
            idleTimer += deltaTime;
            if (idleTimer >= 1.0f)
            {
                TraceLog(
                    LOG_INFO,
                    TextFormat(
                        "jugemu stopped at: (%.3f, %.3f, %.3f)",
                        jugemu.position.x,
                        jugemu.position.y,
                        jugemu.position.z));
                movedSinceLastLog = false;
                idleTimer         = 0.0f;
            }
        }

        //----------------------------------------------------------------------------------
        // Draw
        //----------------------------------------------------------------------------------
        BeginDrawing();
        ClearBackground(BLACK);

        if (PERSPECTIVE_CORRECT())
        {
            capturePerspectiveCorrectToTexture(
                &mainObserver, &mainModel, MODEL_POS, meshRotationRadians, MODEL_SCALE, tex, COLOR_MODE());
            ClearBackground(BLACK);
        }

        BeginMode3D(jugemu);
        drawObservedAxes(&mainObserver);

        if ((gFlags & FLAG_NDC_OVERLAY) ? 1 : 0)
        {
            rlDisableBackfaceCulling();

            mapFrustumToNdcCube(
                &mainObserver,
                screenWidth,
                screenHeight,
                nearClipPlane,
                farClipPlane,
                &mainModel,
                &ndcModel,
                MODEL_POS,
                MODEL_SCALE,
                meshRotationRadians);
            ndcModel.materials[0].maps[MATERIAL_MAP_ALBEDO].texture.id = tex.id;
            ndcModel.materials[0].maps[MATERIAL_MAP_ALBEDO].texture.id = COLOR_MODE() ? 0 : tex.id;
            DrawModelEx(ndcModel, MODEL_POS, (Vector3){0, 1, 0}, RAD2DEG * meshRotationRadians, MODEL_SCALE, WHITE);
            rlSetLineWidth(2.0f);
            ndcModel.materials[0].maps[MATERIAL_MAP_ALBEDO].texture.id = 0;
            DrawModelWiresEx(ndcModel, MODEL_POS, (Vector3){0, 1, 0}, RAD2DEG * meshRotationRadians, MODEL_SCALE, BLUE);
            rlSetPointSize(8.0f);
            DrawModelPointsEx(
                ndcModel, MODEL_POS, (Vector3){0, 1, 0}, RAD2DEG * meshRotationRadians, MODEL_SCALE, GREEN);

            Topology topoNdc;
            topologyBuild(&topoNdc, &ndcModel.meshes[0]);
            topologyFrontTriangles(&topoNdc, meshRotationRadians, &mainObserver);
            topologySilhouetteTriangles(&topoNdc);

            drawNearPlaneIntersectionalDiskMesh(
                &mainObserver,
                nearClipPlane,
                &nearPlaneIntersectionalDiskModel,
                MODEL_POS,
                MODEL_SCALE,
                meshRotationRadians,
                &topoNdc,
                true);

            if (PERSPECTIVE_CORRECT() && gPerspectiveCorrectTextureReady)
            {
                ensureNearPlaneQuadBuilt();
                updateNearPlaneQuadGeometry(&mainObserver, nearClipPlane, screenWidth, screenHeight);
                gNearQuadModel.materials[0].maps[MATERIAL_MAP_ALBEDO].texture = gPerspectiveCorrectTexture;
                rlEnableColorBlend();
                DrawModel(gNearQuadModel, (Vector3){0, 0, 0}, 1.0f, WHITE);
            }
            else
            {
                if (COLOR_MODE())
                {
                    perspectiveIncorrectColorDidactic(
                        screenWidth,
                        screenHeight,
                        &mainObserver,
                        nearClipPlane,
                        &ndcModel.meshes[0],
                        MODEL_POS,
                        MODEL_SCALE,
                        meshRotationRadians);
                }
                else
                {
                    perspectiveIncorrectProjectionDidactic(
                        screenWidth,
                        screenHeight,
                        &mainObserver,
                        nearClipPlane,
                        &ndcModel.meshes[0],
                        MODEL_POS,
                        MODEL_SCALE,
                        meshRotationRadians,
                        tex);
                }
            }

            freeTopology(&topoNdc);
        }
        else
        {
            mainModel.materials[0].maps[MATERIAL_MAP_ALBEDO].texture.id = COLOR_MODE() ? 0 : tex.id;
            DrawModelEx(mainModel, MODEL_POS, (Vector3){0, 1, 0}, RAD2DEG * meshRotationRadians, MODEL_SCALE, WHITE);
            rlSetLineWidth(2.0f);
            mainModel.materials[0].maps[MATERIAL_MAP_ALBEDO].texture.id = 0;
            DrawModelWiresEx(
                mainModel, MODEL_POS, (Vector3){0, 1, 0}, RAD2DEG * meshRotationRadians, MODEL_SCALE, BLUE);
            rlSetPointSize(8.0f);
            DrawModelPointsEx(
                mainModel, MODEL_POS, (Vector3){0, 1, 0}, RAD2DEG * meshRotationRadians, MODEL_SCALE, GREEN);

            drawFrustum(&mainObserver, ANISOTROPIC() ? aspect : 1.0f, nearClipPlane, farClipPlane);

            Topology topo;
            topologyBuild(&topo, &mainModel.meshes[0]);
            topologyFrontTriangles(&topo, meshRotationRadians, &mainObserver);
            topologySilhouetteTriangles(&topo);

            drawNearPlaneIntersectionalDiskMesh(
                &mainObserver,
                nearClipPlane,
                &nearPlaneIntersectionalDiskModel,
                MODEL_POS,
                MODEL_SCALE,
                meshRotationRadians,
                &topo,
                false);

            if (PERSPECTIVE_CORRECT() && gPerspectiveCorrectTextureReady)
            {
                ensureNearPlaneQuadBuilt();
                updateNearPlaneQuadGeometry(&mainObserver, nearClipPlane, screenWidth, screenHeight);
                gNearQuadModel.materials[0].maps[MATERIAL_MAP_ALBEDO].texture = gPerspectiveCorrectTexture;
                rlEnableColorBlend();
                DrawModel(gNearQuadModel, (Vector3){0, 0, 0}, 1.0f, WHITE);
            }
            else
            {
                if (COLOR_MODE())
                {
                    perspectiveIncorrectColorDidactic(
                        screenWidth,
                        screenHeight,
                        &mainObserver,
                        nearClipPlane,
                        &mainModel.meshes[0],
                        MODEL_POS,
                        MODEL_SCALE,
                        meshRotationRadians);
                }
                else
                {
                    perspectiveIncorrectProjectionDidactic(
                        screenWidth,
                        screenHeight,
                        &mainObserver,
                        nearClipPlane,
                        &mainModel.meshes[0],
                        MODEL_POS,
                        MODEL_SCALE,
                        meshRotationRadians,
                        tex);
                }
            }

            freeTopology(&topo);
        }

        // Optional software raster
        // int rasterStep = 4;
        // drawNearPlaneSoftwareRaster(&mainObserver, screenWidth, screenHeight, nearClipPlane,
        //                             &mainModel.meshes[0], MODEL_POS, MODEL_SCALE, meshRotationRadians, rasterStep);

        EndMode3D();

        DrawText(
            "E: NDC overlay  ·  Q: isotropy  ·  P: perspective-correctness  ·  F: freeze  ·  Arrows: orbit  ·  W/S: zoom  ·  A/D: roll  ·  Space: reset",
            12,
            14,
            10,
            RAYWHITE);

        DrawText(
            TextFormat(
                "NDC:%u  Aspect:%s  Proj:%s  Freeze:%u",
                (gFlags & FLAG_NDC_OVERLAY) ? 1 : 0,
                ANISOTROPIC() ? "ANISO" : "ISO",
                PERSPECTIVE_CORRECT() ? "CORRECT" : "INCORRECT",
                (gFlags & FLAG_PAUSE) ? 1 : 0),
            12,
            28,
            10,
            RAYWHITE);

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

//----------------------------------------------------------------------------------
// Small math helpers
//----------------------------------------------------------------------------------
static inline Vector3 observedLineOfSight(Camera3D *observer)
{
    Vector3 lineOfSight = Vector3Subtract(observer->target, observer->position);
    return Vector3Normalize(lineOfSight);
}

static inline void observedBasis(
    Camera3D *observer, Vector3 *observedLineOfSightOut, Vector3 *observedRightOut, Vector3 *observedUpOut)
{
    Vector3 los             = observedLineOfSight(observer);
    Vector3 right           = Vector3Normalize(Vector3CrossProduct(los, observer->up));
    Vector3 up              = Vector3Normalize(Vector3CrossProduct(right, los));
    *observedLineOfSightOut = los;
    *observedRightOut       = right;
    *observedUpOut          = up;
}

static inline Vector3 rotatePointAboutAxis(Vector3 point, Vector3 axisStart, Vector3 axisEnd, float angleRadians)
{
    Vector3 axisDir       = Vector3Normalize(Vector3Subtract(axisEnd, axisStart));
    Vector3 localFromAxis = Vector3Subtract(point, axisStart);
    Vector3 rotatedLocal  = Vector3RotateByAxisAngle(localFromAxis, axisDir, angleRadians);
    return Vector3Add(axisStart, rotatedLocal);
}

static inline void rotateVerticesInPlaneSlice(Vector3 *vertices, int vertexCount, float rotationRadians)
{
    float s = sinf(rotationRadians), c = cosf(rotationRadians);
    for (int vertexIndex = 0; vertexIndex < vertexCount; vertexIndex++)
    {
        float x0 = vertices[vertexIndex].x, z0 = vertices[vertexIndex].z;
        vertices[vertexIndex].x = c * x0 + s * z0;
        vertices[vertexIndex].z = -s * x0 + c * z0;
    }
}

static inline Vector3
    applyModelTranslateRotateScale(Vector3 modelCoord, Vector3 modelPosition, Vector3 modelScale, float rotationRadians)
{
    Matrix M = MatrixMultiply(
        MatrixMultiply(MatrixScale(modelScale.x, modelScale.y, modelScale.z), MatrixRotateY(rotationRadians)),
        MatrixTranslate(modelPosition.x, modelPosition.y, modelPosition.z));
    return Vector3Transform(modelCoord, M);
}

static inline Vector3 applyInverseModelTranslateRotateScale(
    Vector3 worldCoord, Vector3 modelPosition, Vector3 modelScale, float rotationRadians)
{
    Matrix M = MatrixMultiply(
        MatrixMultiply(MatrixScale(modelScale.x, modelScale.y, modelScale.z), MatrixRotateY(rotationRadians)),
        MatrixTranslate(modelPosition.x, modelPosition.y, modelPosition.z));
    Matrix invM = MatrixInvert(M);
    return Vector3Transform(worldCoord, invM);
}

static inline Vector3 nearPlaneIntersection(Camera3D *observer, float nearClipPlane, Vector3 worldCoord)
{
    Vector3 los        = observedLineOfSight(observer);
    Vector3 rayToWorld = Vector3Subtract(worldCoord, observer->position);
    float signedDepth  = Vector3DotProduct(rayToWorld, los); // positive forward depth
    float t            = nearClipPlane / (signedDepth + 1e-9f);
    return Vector3Add(observer->position, Vector3Scale(rayToWorld, t));
}

static inline Vector3 triangleNormal(Vector3 a, Vector3 b, Vector3 c)
{
    return Vector3Normalize(Vector3CrossProduct(Vector3Subtract(b, a), Vector3Subtract(c, a)));
}

// FLAG STUIFF:

// Build a 2-tri quad with UVs flipped vertically to sample a RenderTexture correctly.
static void ensureNearPlaneQuadBuilt(void)
{
    if (gNearQuadBuilt)
        return;

    gNearQuad               = (Mesh){0};
    gNearQuad.vertexCount   = 4;
    gNearQuad.triangleCount = 2;
    gNearQuad.vertices      = (float *)RL_CALLOC(4 * 3, sizeof(float));
    gNearQuad.texcoords     = (float *)RL_CALLOC(4 * 2, sizeof(float));
    gNearQuad.indices       = (unsigned short *)RL_CALLOC(6, sizeof(unsigned short));

    // Standard UVs for CPU screenshot textures (no flip)
    float uvs[8] = {
        0.0f,
        0.0f, // 0 = TL
        1.0f,
        0.0f, // 1 = TR
        1.0f,
        1.0f, // 2 = BR
        0.0f,
        1.0f // 3 = BL
    };
    memcpy(gNearQuad.texcoords, uvs, sizeof(uvs));

    unsigned short idx[6] = {0, 2, 1, 0, 3, 2};
    memcpy(gNearQuad.indices, idx, sizeof(idx));

    UploadMesh(&gNearQuad, false);
    gNearQuadModel = LoadModelFromMesh(gNearQuad);
    gNearQuadBuilt = true;
}

static void updateNearPlaneQuadGeometry(Camera3D *observer, float nearClipPlane, int screenWidth, int screenHeight)
{
    Vector3 los, right, up;
    observedBasis(observer, &los, &right, &up);

    float aspect   = (float)screenWidth / (float)screenHeight;
    float halfFovy = DEG2RAD * observer->fovy * 0.5f;
    float halfH    = nearClipPlane * tanf(halfFovy);
    float halfW    = ANISOTROPIC() ? (halfH * aspect) : halfH;

    Vector3 centerNear = Vector3Add(observer->position, Vector3Scale(los, nearClipPlane));

    Vector3 TL = Vector3Add(centerNear, Vector3Add(Vector3Scale(up, +halfH), Vector3Scale(right, -halfW)));
    Vector3 TR = Vector3Add(centerNear, Vector3Add(Vector3Scale(up, +halfH), Vector3Scale(right, +halfW)));
    Vector3 BR = Vector3Add(centerNear, Vector3Add(Vector3Scale(up, -halfH), Vector3Scale(right, +halfW)));
    Vector3 BL = Vector3Add(centerNear, Vector3Add(Vector3Scale(up, -halfH), Vector3Scale(right, -halfW)));

    float *v = gNearQuad.vertices;
    v[0]     = TL.x;
    v[1]     = TL.y;
    v[2]     = TL.z;
    v[3]     = TR.x;
    v[4]     = TR.y;
    v[5]     = TR.z;
    v[6]     = BR.x;
    v[7]     = BR.y;
    v[8]     = BR.z;
    v[9]     = BL.x;
    v[10]    = BL.y;
    v[11]    = BL.z;

    UploadMesh(&gNearQuad, false);
}

static void composeAlphaFromMask(Image *rgba, const Image *mask, unsigned char threshold)
{
    // Preconditions: same size
    if (rgba->width != mask->width || rgba->height != mask->height)
        return;

    ImageFormat(rgba, PIXELFORMAT_UNCOMPRESSED_R8G8B8A8);
    Image maskCopy = ImageCopy(*mask);
    ImageFormat(&maskCopy, PIXELFORMAT_UNCOMPRESSED_R8G8B8A8);

    unsigned char *dst = (unsigned char *)rgba->data;
    unsigned char *msk = (unsigned char *)maskCopy.data;
    int n              = rgba->width * rgba->height;

    for (int i = 0; i < n; ++i)
    {
        unsigned int v  = (unsigned int)msk[4 * i + 0] + (unsigned int)msk[4 * i + 1] + (unsigned int)msk[4 * i + 2];
        unsigned char a = (unsigned char)(v / 3);
        dst[4 * i + 3]  = (a > threshold) ? a : 0;
    }

    UnloadImage(maskCopy);
}

static void capturePerspectiveCorrectToTexture(
    Camera3D *mainObserver,
    Model *model,
    Vector3 modelPos,
    float rotRadians,
    Vector3 modelScale,
    Texture2D tex,
    bool usePVC)
{
    // ---------- PASS 1: normal color ----------
    ClearBackground(BLACK);
    BeginMode3D(*mainObserver);
    Texture2D prevTex                                     = model->materials[0].maps[MATERIAL_MAP_ALBEDO].texture;
    model->materials[0].maps[MATERIAL_MAP_ALBEDO].texture = usePVC ? (Texture2D){0} : tex;
    DrawModelEx(*model, modelPos, (Vector3){0, 1, 0}, RAD2DEG * rotRadians, modelScale, WHITE);
    model->materials[0].maps[MATERIAL_MAP_ALBEDO].texture = prevTex;
    EndMode3D();

    Image colorShot = LoadImageFromScreen();
    ImageFormat(&colorShot, PIXELFORMAT_UNCOMPRESSED_R8G8B8A8);

    // ---------- PASS 2: coverage mask (solid white mesh on black) ----------
    ClearBackground(BLACK);
    BeginMode3D(*mainObserver);
    Texture2D saveTex = model->materials[0].maps[MATERIAL_MAP_ALBEDO].texture;
    Color saveCol     = model->materials[0].maps[MATERIAL_MAP_ALBEDO].color;

    model->materials[0].maps[MATERIAL_MAP_ALBEDO].texture = (Texture2D){0};
    model->materials[0].maps[MATERIAL_MAP_ALBEDO].color   = WHITE;

    DrawModelEx(*model, modelPos, (Vector3){0, 1, 0}, RAD2DEG * rotRadians, modelScale, WHITE);

    model->materials[0].maps[MATERIAL_MAP_ALBEDO].texture = saveTex;
    model->materials[0].maps[MATERIAL_MAP_ALBEDO].color   = saveCol;
    EndMode3D();

    Image maskShot = LoadImageFromScreen();
    composeAlphaFromMask(&colorShot, &maskShot, 1);

    if (gPerspectiveCorrectTextureReady)
    {
        UpdateTexture(gPerspectiveCorrectTexture, colorShot.data);
    }
    else
    {
        gPerspectiveCorrectTexture      = LoadTextureFromImage(colorShot);
        gPerspectiveCorrectTextureReady = true;
    }

    UnloadImage(maskShot);
    UnloadImage(colorShot);
}

//----------------------------------------------------------------------------------
// Drawing helpers
//----------------------------------------------------------------------------------
static void drawObservedAxes(Camera3D *observer)
{
    Vector3 los, right, up;
    observedBasis(observer, &los, &right, &up);
    rlSetLineWidth(1.5f);
    DrawLine3D(observer->position, Vector3Add(observer->position, right), PURPLE);
    DrawLine3D(observer->position, Vector3Add(observer->position, up), DARKGREEN);
    DrawLine3D(observer->position, Vector3Add(observer->position, los), SKYBLUE);
}

static void drawFrustum(Camera3D *observer, float aspect, float nearClipPlane, float farClipPlane)
{
    Vector3 los, right, up;
    observedBasis(observer, &los, &right, &up);

    Vector3 centerNear = Vector3Add(observer->position, Vector3Scale(los, nearClipPlane));
    Vector3 centerFar  = Vector3Add(observer->position, Vector3Scale(los, farClipPlane));

    float halfFovy       = DEG2RAD * observer->fovy * 0.5f;
    float halfHeightNear = nearClipPlane * tanf(halfFovy);
    float halfWidthNear  = halfHeightNear * aspect;
    float halfHeightFar  = farClipPlane * tanf(halfFovy);
    float halfWidthFar   = halfHeightFar * aspect;

    Vector3 nearTL =
        Vector3Add(Vector3Add(centerNear, Vector3Scale(up, halfHeightNear)), Vector3Scale(right, -halfWidthNear));
    Vector3 nearTR =
        Vector3Add(Vector3Add(centerNear, Vector3Scale(up, halfHeightNear)), Vector3Scale(right, halfWidthNear));
    Vector3 nearBR =
        Vector3Add(Vector3Add(centerNear, Vector3Scale(up, -halfHeightNear)), Vector3Scale(right, halfWidthNear));
    Vector3 nearBL =
        Vector3Add(Vector3Add(centerNear, Vector3Scale(up, -halfHeightNear)), Vector3Scale(right, -halfWidthNear));

    Vector3 farTL =
        Vector3Add(Vector3Add(centerFar, Vector3Scale(up, halfHeightFar)), Vector3Scale(right, -halfWidthFar));
    Vector3 farTR =
        Vector3Add(Vector3Add(centerFar, Vector3Scale(up, halfHeightFar)), Vector3Scale(right, halfWidthFar));
    Vector3 farBR =
        Vector3Add(Vector3Add(centerFar, Vector3Scale(up, -halfHeightFar)), Vector3Scale(right, halfWidthFar));
    Vector3 farBL =
        Vector3Add(Vector3Add(centerFar, Vector3Scale(up, -halfHeightFar)), Vector3Scale(right, -halfWidthFar));

    rlSetLineWidth(1.0f);
    DrawLine3D(nearTL, nearTR, SKYBLUE);
    DrawLine3D(nearTR, nearBR, SKYBLUE);
    DrawLine3D(nearBR, nearBL, SKYBLUE);
    DrawLine3D(nearBL, nearTL, SKYBLUE);

    DrawLine3D(farTL, farTR, GRAY);
    DrawLine3D(farTR, farBR, GRAY);
    DrawLine3D(farBR, farBL, GRAY);
    DrawLine3D(farBL, farTL, GRAY);

    DrawLine3D(nearTL, farTL, DARKBLUE);
    DrawLine3D(nearTR, farTR, DARKBLUE);
    DrawLine3D(nearBR, farBR, DARKBLUE);
    DrawLine3D(nearBL, farBL, DARKBLUE);
}

static void drawNdcCube(
    Camera3D *observer,
    Vector3 centerNear,
    float halfWidthNear,
    float halfHeightNear,
    float halfDepthNdcCube,
    Color edgeColor)
{
    Vector3 los, right, up;
    observedBasis(observer, &los, &right, &up);

    Vector3 corners[8];
    int writeIndex = 0;
    for (int zSign = -1; zSign <= 1; zSign += 2)
        for (int ySign = -1; ySign <= 1; ySign += 2)
            for (int xSign = -1; xSign <= 1; xSign += 2)
            {
                corners[writeIndex++] = Vector3Add(
                    centerNear,
                    Vector3Add(
                        Vector3Scale(right, xSign * halfWidthNear),
                        Vector3Add(
                            Vector3Scale(up, ySign * halfHeightNear), Vector3Scale(los, zSign * halfDepthNdcCube))));
            }

    Vector3 nearTL = corners[3], nearTR = corners[2], nearBR = corners[0], nearBL = corners[1];
    Vector3 farTL = corners[7], farTR = corners[6], farBR = corners[4], farBL = corners[5];

    rlSetLineWidth(1.0f);
    DrawLine3D(nearTL, nearTR, edgeColor);
    DrawLine3D(nearTR, nearBR, edgeColor);
    DrawLine3D(nearBR, nearBL, edgeColor);
    DrawLine3D(nearBL, nearTL, edgeColor);
    DrawLine3D(farTL, farTR, edgeColor);
    DrawLine3D(farTR, farBR, edgeColor);
    DrawLine3D(farBR, farBL, edgeColor);
    DrawLine3D(farBL, farTL, edgeColor);
    DrawLine3D(nearTL, farTL, edgeColor);
    DrawLine3D(nearTR, farTR, edgeColor);
    DrawLine3D(nearBR, farBR, edgeColor);
    DrawLine3D(nearBL, farBL, edgeColor);
}

//----------------------------------------------------------------------------------
// World -> NDC mapping
//----------------------------------------------------------------------------------
static Vector3
    worldCoordToNdcCoord(float aspect, Camera3D *observer, float nearClipPlane, float farClipPlane, Vector3 worldCoord)
{
    Vector3 los, right, up;
    observedBasis(observer, &los, &right, &up);

    float halfFovy       = DEG2RAD * observer->fovy * 0.5f;
    float halfHeightNear = nearClipPlane * tanf(halfFovy);
    float halfWidthNear  = halfHeightNear * aspect;

    Vector3 intersectionCoord = nearPlaneIntersection(observer, nearClipPlane, worldCoord);
    Vector3 centerNear        = Vector3Add(observer->position, Vector3Scale(los, nearClipPlane));
    Vector3 clipPlaneVector   = Vector3Subtract(intersectionCoord, centerNear);

    float xNdc = Vector3DotProduct(clipPlaneVector, right) / (halfWidthNear + 1e-9f);
    float yNdc = Vector3DotProduct(clipPlaneVector, up) / (halfHeightNear + 1e-9f);

    float signedDepth = Vector3DotProduct(Vector3Subtract(worldCoord, observer->position), los);
    float zNdc = ((farClipPlane + nearClipPlane) - (2.0f * farClipPlane * nearClipPlane) / (signedDepth + 1e-9f)) /
                 (farClipPlane - nearClipPlane);

    return (Vector3){xNdc, yNdc, zNdc};
}

static Vector3 scaleNdcCoordByNearClipPlane(
    Camera3D *observer,
    Vector3 centerNear,
    float halfWidthNear,
    float halfHeightNear,
    float halfDepthNdcCube,
    Vector3 ndcCoord)
{
    Vector3 los, right, up;
    observedBasis(observer, &los, &right, &up);
    return Vector3Add(
        centerNear,
        Vector3Add(
            Vector3Scale(right, ndcCoord.x * halfWidthNear),
            Vector3Add(
                Vector3Scale(up, ndcCoord.y * halfHeightNear), Vector3Scale(los, ndcCoord.z * halfDepthNdcCube))));
}

static void updateWorldToNdcMappedMesh(
    Mesh *ndcMesh,
    Mesh *worldMesh,
    int screenWidth,
    int screenHeight,
    Camera3D *observer,
    float nearClipPlane,
    float farClipPlane,
    Vector3 modelPosition,
    Vector3 modelScale,
    float meshRotationRadians)
{
    float aspect = (float)screenWidth / (float)screenHeight;
    Vector3 los  = observedLineOfSight(observer);

    float halfFovy       = DEG2RAD * observer->fovy * 0.5f;
    float halfHeightNear = nearClipPlane * tanf(halfFovy);

    float halfWidthNear = ANISOTROPIC() ? (halfHeightNear * aspect) : halfHeightNear;
    float halfDepthNdc  = ANISOTROPIC() ? (0.5f * (farClipPlane - nearClipPlane)) : halfHeightNear;

    Vector3 centerNear     = Vector3Add(observer->position, Vector3Scale(los, nearClipPlane));
    Vector3 ndcCubeCenter  = Vector3Add(centerNear, Vector3Scale(los, halfDepthNdc));
    float aspectForMapping = ANISOTROPIC() ? aspect : 1.0f;
    if (worldMesh->indices && ndcMesh->indices)
    {
        memcpy(ndcMesh->indices, worldMesh->indices, sizeof(unsigned short) * 3 * worldMesh->triangleCount);
        ndcMesh->triangleCount = worldMesh->triangleCount;
    }

    for (int vertexIndex = 0; vertexIndex < worldMesh->vertexCount; vertexIndex++)
    {
        Vector3 objectVertex = {
            worldMesh->vertices[3 * vertexIndex + 0],
            worldMesh->vertices[3 * vertexIndex + 1],
            worldMesh->vertices[3 * vertexIndex + 2]};
        Vector3 worldVertex =
            applyModelTranslateRotateScale(objectVertex, modelPosition, modelScale, meshRotationRadians);
        Vector3 ndcCoord = worldCoordToNdcCoord(aspectForMapping, observer, nearClipPlane, farClipPlane, worldVertex);
        Vector3 scaledNdcCoord = scaleNdcCoordByNearClipPlane(
            observer, ndcCubeCenter, halfWidthNear, halfHeightNear, halfDepthNdc, ndcCoord);
        Vector3 mappedObjectCoord =
            applyInverseModelTranslateRotateScale(scaledNdcCoord, modelPosition, modelScale, meshRotationRadians);

        ndcMesh->vertices[3 * vertexIndex + 0] = mappedObjectCoord.x;
        ndcMesh->vertices[3 * vertexIndex + 1] = mappedObjectCoord.y;
        ndcMesh->vertices[3 * vertexIndex + 2] = mappedObjectCoord.z;
    }

    UploadMesh(ndcMesh, false);
}

static void mapFrustumToNdcCube(
    Camera3D *observer,
    int screenWidth,
    int screenHeight,
    float nearClipPlane,
    float farClipPlane,
    Model *worldModel,
    Model *ndcModel,
    Vector3 modelPosition,
    Vector3 modelScale,
    float meshRotationRadians)

{
    float aspect         = (float)screenWidth / (float)screenHeight;
    Vector3 los          = observedLineOfSight(observer);
    float halfFovy       = DEG2RAD * observer->fovy * 0.5f;
    float halfHeightNear = nearClipPlane * tanf(halfFovy);

    float halfWidthNear = ANISOTROPIC() ? (halfHeightNear * aspect) : halfHeightNear;
    float halfDepthNdc  = ANISOTROPIC() ? (0.5f * (farClipPlane - nearClipPlane)) : halfHeightNear;

    Vector3 centerNear    = Vector3Add(observer->position, Vector3Scale(los, nearClipPlane));
    Vector3 ndcCubeCenter = Vector3Add(centerNear, Vector3Scale(los, halfDepthNdc));

    drawNdcCube(observer, ndcCubeCenter, halfWidthNear, halfHeightNear, halfDepthNdc, SKYBLUE);

    updateWorldToNdcMappedMesh(
        &ndcModel->meshes[0],
        &worldModel->meshes[0],
        screenWidth,
        screenHeight,
        observer,
        nearClipPlane,
        farClipPlane,
        modelPosition,
        modelScale,
        meshRotationRadians);
}

static void freeTopology(Topology *topology)
{
    if (!topology)
        return;
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
    if (topology->triangleCount <= 0)
        return;

    topology->triangles = (unsigned short *)malloc(sizeof(unsigned short) * 3 * topology->triangleCount);
    if (mesh->indices)
        memcpy(topology->triangles, mesh->indices, sizeof(unsigned short) * 3 * topology->triangleCount);
    else
        for (int linearIndex = 0; linearIndex < 3 * topology->triangleCount; linearIndex++)
            topology->triangles[linearIndex] = (unsigned short)linearIndex;

    topology->verticesPerTriangle = (Vector3 *)malloc(sizeof(Vector3) * 3 * topology->triangleCount);
    for (int triangleIndex = 0; triangleIndex < topology->triangleCount; triangleIndex++)
    {
        unsigned short indexA = topology->triangles[3 * triangleIndex + 0];
        unsigned short indexB = topology->triangles[3 * triangleIndex + 1];
        unsigned short indexC = topology->triangles[3 * triangleIndex + 2];
        topology->verticesPerTriangle[3 * triangleIndex + 0] =
            (Vector3){mesh->vertices[3 * indexA + 0], mesh->vertices[3 * indexA + 1], mesh->vertices[3 * indexA + 2]};
        topology->verticesPerTriangle[3 * triangleIndex + 1] =
            (Vector3){mesh->vertices[3 * indexB + 0], mesh->vertices[3 * indexB + 1], mesh->vertices[3 * indexB + 2]};
        topology->verticesPerTriangle[3 * triangleIndex + 2] =
            (Vector3){mesh->vertices[3 * indexC + 0], mesh->vertices[3 * indexC + 1], mesh->vertices[3 * indexC + 2]};
    }

    topology->weldedVertexId = (int *)malloc(sizeof(int) * mesh->vertexCount);
    for (int vertexIndex = 0; vertexIndex < mesh->vertexCount; vertexIndex++)
    {
        int quantX      = (int)lroundf(mesh->vertices[3 * vertexIndex + 0] / 1e-5f);
        int quantY      = (int)lroundf(mesh->vertices[3 * vertexIndex + 1] / 1e-5f);
        int quantZ      = (int)lroundf(mesh->vertices[3 * vertexIndex + 2] / 1e-5f);
        int foundWelded = -1;
        for (int previous = 0; previous < vertexIndex; previous++)
        {
            int quantXPrev = (int)lroundf(mesh->vertices[3 * previous + 0] / 1e-5f);
            int quantYPrev = (int)lroundf(mesh->vertices[3 * previous + 1] / 1e-5f);
            int quantZPrev = (int)lroundf(mesh->vertices[3 * previous + 2] / 1e-5f);
            if (quantX == quantXPrev && quantY == quantYPrev && quantZ == quantZPrev)
            {
                foundWelded = topology->weldedVertexId[previous];
                break;
            }
        }
        topology->weldedVertexId[vertexIndex] = (foundWelded >= 0) ? foundWelded : vertexIndex;
    }

    topology->weldedVerticesPerTriangle = (WeldedVertex *)malloc(sizeof(WeldedVertex) * 3 * topology->triangleCount);
    for (int triangleIndex = 0; triangleIndex < topology->triangleCount; triangleIndex++)
    {
        unsigned short indexA                                         = topology->triangles[3 * triangleIndex + 0];
        unsigned short indexB                                         = topology->triangles[3 * triangleIndex + 1];
        unsigned short indexC                                         = topology->triangles[3 * triangleIndex + 2];
        topology->weldedVerticesPerTriangle[3 * triangleIndex + 0].id = topology->weldedVertexId[indexA];
        topology->weldedVerticesPerTriangle[3 * triangleIndex + 1].id = topology->weldedVertexId[indexB];
        topology->weldedVerticesPerTriangle[3 * triangleIndex + 2].id = topology->weldedVertexId[indexC];
    }

    topology->neighborsPerTriangle = (int (*)[3])malloc(sizeof(int) * 3 * topology->triangleCount);
    for (int triangleIndex = 0; triangleIndex < topology->triangleCount; triangleIndex++)
        topology->neighborsPerTriangle[triangleIndex][0]     = topology->neighborsPerTriangle[triangleIndex][1] =
            topology->neighborsPerTriangle[triangleIndex][2] = -1;

    for (int triangleIndex = 0; triangleIndex < topology->triangleCount; triangleIndex++)
    {
        int edge0[2] = {
            topology->weldedVerticesPerTriangle[3 * triangleIndex + 0].id,
            topology->weldedVerticesPerTriangle[3 * triangleIndex + 1].id};
        int edge1[2] = {
            topology->weldedVerticesPerTriangle[3 * triangleIndex + 1].id,
            topology->weldedVerticesPerTriangle[3 * triangleIndex + 2].id};
        int edge2[2] = {
            topology->weldedVerticesPerTriangle[3 * triangleIndex + 2].id,
            topology->weldedVerticesPerTriangle[3 * triangleIndex + 0].id};

        for (int otherTriangleIndex = triangleIndex + 1; otherTriangleIndex < topology->triangleCount;
             otherTriangleIndex++)
        {
            int otherEdge0[2] = {
                topology->weldedVerticesPerTriangle[3 * otherTriangleIndex + 0].id,
                topology->weldedVerticesPerTriangle[3 * otherTriangleIndex + 1].id};
            int otherEdge1[2] = {
                topology->weldedVerticesPerTriangle[3 * otherTriangleIndex + 1].id,
                topology->weldedVerticesPerTriangle[3 * otherTriangleIndex + 2].id};
            int otherEdge2[2] = {
                topology->weldedVerticesPerTriangle[3 * otherTriangleIndex + 2].id,
                topology->weldedVerticesPerTriangle[3 * otherTriangleIndex + 0].id};

            int weldedEdgePairsThis[3][2]  = {{edge0[0], edge0[1]}, {edge1[0], edge1[1]}, {edge2[0], edge2[1]}};
            int weldedEdgePairsOther[3][2] = {
                {otherEdge0[0], otherEdge0[1]}, {otherEdge1[0], otherEdge1[1]}, {otherEdge2[0], otherEdge2[1]}};

            for (int edgeIndexThis = 0; edgeIndexThis < 3; edgeIndexThis++)
                for (int edgeIndexOther = 0; edgeIndexOther < 3; edgeIndexOther++)
                {
                    int edge0VMin = weldedEdgePairsThis[edgeIndexThis][0];
                    int edge0VMax = weldedEdgePairsThis[edgeIndexThis][1];
                    if (edge0VMin > edge0VMax)
                    {
                        int swapVertex = edge0VMin;
                        edge0VMin      = edge0VMax;
                        edge0VMax      = swapVertex;
                    }

                    int edge1VMin = weldedEdgePairsOther[edgeIndexOther][0];
                    int edge1VMax = weldedEdgePairsOther[edgeIndexOther][1];
                    if (edge1VMin > edge1VMax)
                    {
                        int swapVertex = edge1VMin;
                        edge1VMin      = edge1VMax;
                        edge1VMax      = swapVertex;
                    }

                    if (edge0VMin == edge1VMin && edge0VMax == edge1VMax)
                    {
                        topology->neighborsPerTriangle[triangleIndex][edgeIndexThis]       = otherTriangleIndex;
                        topology->neighborsPerTriangle[otherTriangleIndex][edgeIndexOther] = triangleIndex;
                    }
                }
        }
    }

    topology->frontTrianglesFlag      = (unsigned char *)calloc(topology->triangleCount, 1);
    topology->silhouetteTrianglesFlag = (unsigned char *)calloc(topology->triangleCount, 1);
}

static void topologyFrontTriangles(Topology *topology, float rotationRadians, Camera3D *observer)
{
    Vector3 los = observedLineOfSight(observer);
    for (int triangleIndex = 0; triangleIndex < topology->triangleCount; triangleIndex++)
    {
        Vector3 triVerts[3] = {
            topology->verticesPerTriangle[3 * triangleIndex + 0],
            topology->verticesPerTriangle[3 * triangleIndex + 1],
            topology->verticesPerTriangle[3 * triangleIndex + 2]};
        rotateVerticesInPlaneSlice(triVerts, 3, rotationRadians);
        Vector3 normal = triangleNormal(triVerts[0], triVerts[1], triVerts[2]);
        if (Vector3DotProduct(normal, los) <= 0.0f)
            topology->frontTrianglesFlag[triangleIndex] = 1;
    }
}

static void topologySilhouetteTriangles(Topology *topology)
{
    for (int triangleIndex = 0; triangleIndex < topology->triangleCount; triangleIndex++)
    {
        for (int edgeIndex = 0; edgeIndex < 3; edgeIndex++)
        {
            int neighborTriangle = topology->neighborsPerTriangle[triangleIndex][edgeIndex];
            if (neighborTriangle < 0)
                continue;
            if (topology->frontTrianglesFlag[triangleIndex] != topology->frontTrianglesFlag[neighborTriangle])
            {
                topology->silhouetteTrianglesFlag[triangleIndex]    = 1;
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
    Vector3 los, right, up;
    observedBasis(observer, &los, &right, &up);
    Vector3 centerNear  = Vector3Add(observer->position, Vector3Scale(los, nearClipPlane));
    Vector3 toClipPlane = Vector3Subtract(intersectionCoord, centerNear);
    float xComponent    = Vector3DotProduct(toClipPlane, right);
    float yComponent    = Vector3DotProduct(toClipPlane, up);
    return Vector3Add(centerNear, Vector3Add(Vector3Scale(right, xComponent), Vector3Scale(up, -yComponent)));
}

static void drawNearPlaneIntersectionalDiskMesh(
    Camera3D *observer,
    float nearClipPlane,
    Model *nearPlaneIntersectionalDiskModel,
    Vector3 modelPosition,
    Vector3 modelScale,
    float meshRotationRadians,
    Topology *topology,
    bool reflectYAxis)
{
    Mesh *diskMesh                  = &nearPlaneIntersectionalDiskModel->meshes[0];
    int diskMeshCapacityVertexCount = diskMesh->vertexCount;

    int intersectionWriteIndex = 0;
    for (int triangleIndex = 0; triangleIndex < topology->triangleCount; triangleIndex++)
    {
        if (!topology->frontTrianglesFlag[triangleIndex])
            continue;

        Vector3 triVerts[3] = {
            topology->verticesPerTriangle[3 * triangleIndex + 0],
            topology->verticesPerTriangle[3 * triangleIndex + 1],
            topology->verticesPerTriangle[3 * triangleIndex + 2]};
        rotateVerticesInPlaneSlice(triVerts, 3, meshRotationRadians);

        for (int localVertexIndex = 0; localVertexIndex < 3; localVertexIndex++)
        {
            Vector3 scaledVertex = {
                triVerts[localVertexIndex].x * modelScale.x,
                triVerts[localVertexIndex].y * modelScale.y,
                triVerts[localVertexIndex].z * modelScale.z};
            Vector3 worldVertex       = Vector3Add(modelPosition, scaledVertex);
            Vector3 intersectionCoord = nearPlaneIntersection(observer, nearClipPlane, worldVertex);
            if (reflectYAxis)
                intersectionCoord = flipYInNearClipPlane(intersectionCoord, observer, nearClipPlane);

            rlSetLineWidth(1.0f);
            DrawLine3D(worldVertex, intersectionCoord, (Color){255, 0, 0, 80});

            if (intersectionWriteIndex < diskMeshCapacityVertexCount)
            {
                diskMesh->vertices[3 * intersectionWriteIndex + 0] = intersectionCoord.x;
                diskMesh->vertices[3 * intersectionWriteIndex + 1] = intersectionCoord.y;
                diskMesh->vertices[3 * intersectionWriteIndex + 2] = intersectionCoord.z;
                intersectionWriteIndex++;
            }
        }
    }

    if (intersectionWriteIndex > 0)
    {
        diskMesh->vertexCount = intersectionWriteIndex;
        UploadMesh(diskMesh, false);
    }
    rlSetPointSize(6.0f);
    DrawModelPoints(*nearPlaneIntersectionalDiskModel, (Vector3){0, 0, 0}, 1.0f, GREEN);
}

static inline Vector3 remapNearPlaneByAspect(
    Vector3 hit,
    Vector3 centerNear,
    Vector3 right,
    Vector3 up,
    float xScale)
{
    Vector3 d = Vector3Subtract(hit, centerNear);
    float xr  = Vector3DotProduct(d, right);
    float yu  = Vector3DotProduct(d, up);
    return Vector3Add(centerNear, Vector3Add(Vector3Scale(right, xr * xScale), Vector3Scale(up, yu)));
}

static void perspectiveIncorrectColorDidactic(
    int screenWidth,
    int screenHeight,
    Camera3D *observer,
    float nearClipPlane,
    Mesh *mesh,
    Vector3 modelPosition,
    Vector3 modelScale,
    float meshRotationRadians)
{
    if (!mesh->colors)
    {
        TraceLog(LOG_ERROR, "perspectiveIncorrectColorDidactic: mesh->colors is NULL. Fill per-vertex colors first.");
        return;
    }

    float aspect = (float)screenWidth / (float)screenHeight;

    Vector3 los, right, up;
    observedBasis(observer, &los, &right, &up);

    Vector3 centerNear   = Vector3Add(observer->position, Vector3Scale(los, nearClipPlane));
    float halfFovy       = DEG2RAD * observer->fovy * 0.5f;
    float halfHeightNear = nearClipPlane * tanf(halfFovy);
    (void)halfHeightNear;

    float xScale         = ANISOTROPIC() ? 1.0f : (1.0f / aspect);
    const float shiftEps = 0.0001f;

    rlDisableTexture();
    rlEnableColorBlend();
    rlBegin(RL_TRIANGLES);

    for (int tri = 0; tri < mesh->triangleCount; tri++)
    {
        int ia = mesh->indices ? mesh->indices[3 * tri + 0] : 3 * tri + 0;
        int ib = mesh->indices ? mesh->indices[3 * tri + 1] : 3 * tri + 1;
        int ic = mesh->indices ? mesh->indices[3 * tri + 2] : 3 * tri + 2;

        Vector3 A = {mesh->vertices[3 * ia + 0], mesh->vertices[3 * ia + 1], mesh->vertices[3 * ia + 2]};
        Vector3 B = {mesh->vertices[3 * ib + 0], mesh->vertices[3 * ib + 1], mesh->vertices[3 * ib + 2]};
        Vector3 C = {mesh->vertices[3 * ic + 0], mesh->vertices[3 * ic + 1], mesh->vertices[3 * ic + 2]};

        A = applyModelTranslateRotateScale(A, modelPosition, modelScale, meshRotationRadians);
        B = applyModelTranslateRotateScale(B, modelPosition, modelScale, meshRotationRadians);
        C = applyModelTranslateRotateScale(C, modelPosition, modelScale, meshRotationRadians);

        Vector3 hitA = nearPlaneIntersection(observer, nearClipPlane, A);
        Vector3 hitB = nearPlaneIntersection(observer, nearClipPlane, B);
        Vector3 hitC = nearPlaneIntersection(observer, nearClipPlane, C);

        Vector3 planeApos = remapNearPlaneByAspect(hitA, centerNear, right, up, xScale);
        Vector3 planeBpos = remapNearPlaneByAspect(hitB, centerNear, right, up, xScale);
        Vector3 planeCpos = remapNearPlaneByAspect(hitC, centerNear, right, up, xScale);

        Vector3 planeA = Vector3Add(planeApos, Vector3Scale(los, shiftEps));
        Vector3 planeB = Vector3Add(planeBpos, Vector3Scale(los, shiftEps));
        Vector3 planeC = Vector3Add(planeCpos, Vector3Scale(los, shiftEps));

        const unsigned char *cA = &mesh->colors[4 * ia];
        const unsigned char *cB = &mesh->colors[4 * ib];
        const unsigned char *cC = &mesh->colors[4 * ic];

        rlColor4ub(cA[0], cA[1], cA[2], cA[3]);
        rlVertex3f(planeA.x, planeA.y, planeA.z);
        rlColor4ub(cB[0], cB[1], cB[2], cB[3]);
        rlVertex3f(planeB.x, planeB.y, planeB.z);
        rlColor4ub(cC[0], cC[1], cC[2], cC[3]);
        rlVertex3f(planeC.x, planeC.y, planeC.z);
    }

    rlEnd();
}

//----------------------------------------------------------------------------------
// Intentionally incorrect projection (affine UV on near plane)
//----------------------------------------------------------------------------------
static void perspectiveIncorrectProjectionDidactic(
    int screenWidth,
    int screenHeight,
    Camera3D *observer,
    float nearClipPlane,
    Mesh *mesh,
    Vector3 modelPosition,
    Vector3 modelScale,
    float meshRotationRadians,
    Texture2D texture)
{
    float aspect = (float)screenWidth / (float)screenHeight;
    Vector3 los, right, up;
    observedBasis(observer, &los, &right, &up);

    Vector3 centerNear   = Vector3Add(observer->position, Vector3Scale(los, nearClipPlane));
    float halfFovy       = DEG2RAD * observer->fovy * 0.5f;
    float halfHeightNear = nearClipPlane * tanf(halfFovy);
    float halfWidthNear  = ANISOTROPIC() ? (halfHeightNear * aspect) : halfHeightNear;

    float xScale = ANISOTROPIC() ? 1.0f : (1.0f / aspect); // squeeze X in ISO mode

    rlColor4f(1, 1, 1, 1);
    rlEnableTexture(texture.id);
    rlSetTexture(texture.id);
    rlBegin(RL_TRIANGLES);

    for (int triangleIndex = 0; triangleIndex < mesh->triangleCount; triangleIndex++)
    {
        int indexA = mesh->indices ? mesh->indices[3 * triangleIndex + 0] : 3 * triangleIndex + 0;
        int indexB = mesh->indices ? mesh->indices[3 * triangleIndex + 1] : 3 * triangleIndex + 1;
        int indexC = mesh->indices ? mesh->indices[3 * triangleIndex + 2] : 3 * triangleIndex + 2;

        Vector3 vertexA = {
            mesh->vertices[3 * indexA + 0], mesh->vertices[3 * indexA + 1], mesh->vertices[3 * indexA + 2]};
        Vector3 vertexB = {
            mesh->vertices[3 * indexB + 0], mesh->vertices[3 * indexB + 1], mesh->vertices[3 * indexB + 2]};
        Vector3 vertexC = {
            mesh->vertices[3 * indexC + 0], mesh->vertices[3 * indexC + 1], mesh->vertices[3 * indexC + 2]};

        vertexA = applyModelTranslateRotateScale(vertexA, modelPosition, modelScale, meshRotationRadians);
        vertexB = applyModelTranslateRotateScale(vertexB, modelPosition, modelScale, meshRotationRadians);
        vertexC = applyModelTranslateRotateScale(vertexC, modelPosition, modelScale, meshRotationRadians);

        Vector3 hitA = nearPlaneIntersection(observer, nearClipPlane, vertexA);
        Vector3 hitB = nearPlaneIntersection(observer, nearClipPlane, vertexB);
        Vector3 hitC = nearPlaneIntersection(observer, nearClipPlane, vertexC);

        const float shiftEps = 0.0001f;
        Vector3 planeApos    = remapNearPlaneByAspect(hitA, centerNear, right, up, xScale);
        Vector3 planeBpos    = remapNearPlaneByAspect(hitB, centerNear, right, up, xScale);
        Vector3 planeCpos    = remapNearPlaneByAspect(hitC, centerNear, right, up, xScale);

        Vector3 planeA = Vector3Add(planeApos, Vector3Scale(los, shiftEps));
        Vector3 planeB = Vector3Add(planeBpos, Vector3Scale(los, shiftEps));
        Vector3 planeC = Vector3Add(planeCpos, Vector3Scale(los, shiftEps));

        float uA = mesh->texcoords ? mesh->texcoords[2 * indexA + 0] : 0.0f;
        float vA = mesh->texcoords ? mesh->texcoords[2 * indexA + 1] : 0.0f;
        float uB = mesh->texcoords ? mesh->texcoords[2 * indexB + 0] : 0.0f;
        float vB = mesh->texcoords ? mesh->texcoords[2 * indexB + 1] : 0.0f;
        float uC = mesh->texcoords ? mesh->texcoords[2 * indexC + 0] : 0.0f;
        float vC = mesh->texcoords ? mesh->texcoords[2 * indexC + 1] : 0.0f;

        rlTexCoord2f(uA, vA);
        rlVertex3f(planeA.x, planeA.y, planeA.z);
        rlTexCoord2f(uB, vB);
        rlVertex3f(planeB.x, planeB.y, planeB.z);
        rlTexCoord2f(uC, vC);
        rlVertex3f(planeC.x, planeC.y, planeC.z);
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
        mesh->colors = (unsigned char *)RL_CALLOC(mesh->vertexCount * 4, sizeof(unsigned char));
    }
    for (int triangleIndex = 0; triangleIndex < mesh->triangleCount; triangleIndex++)
    {
        int indexA = mesh->indices ? mesh->indices[3 * triangleIndex + 0] : 3 * triangleIndex + 0;
        int indexB = mesh->indices ? mesh->indices[3 * triangleIndex + 1] : 3 * triangleIndex + 1;
        int indexC = mesh->indices ? mesh->indices[3 * triangleIndex + 2] : 3 * triangleIndex + 2;

        mesh->colors[4 * indexA + 0] = 255;
        mesh->colors[4 * indexA + 1] = 0;
        mesh->colors[4 * indexA + 2] = 0;
        mesh->colors[4 * indexA + 3] = 255;
        mesh->colors[4 * indexB + 0] = 0;
        mesh->colors[4 * indexB + 1] = 255;
        mesh->colors[4 * indexB + 2] = 0;
        mesh->colors[4 * indexB + 3] = 255;
        mesh->colors[4 * indexC + 0] = 0;
        mesh->colors[4 * indexC + 1] = 0;
        mesh->colors[4 * indexC + 2] = 255;
        mesh->colors[4 * indexC + 3] = 255;
    }

    if (mesh->texcoords == NULL)
    {
        mesh->texcoords = (float *)RL_CALLOC(mesh->vertexCount * 2, sizeof(float));
    }
    for (int triangleIndex = 0; triangleIndex < mesh->triangleCount; triangleIndex++)
    {
        int indexA = mesh->indices ? mesh->indices[3 * triangleIndex + 0] : 3 * triangleIndex + 0;
        int indexB = mesh->indices ? mesh->indices[3 * triangleIndex + 1] : 3 * triangleIndex + 1;
        int indexC = mesh->indices ? mesh->indices[3 * triangleIndex + 2] : 3 * triangleIndex + 2;

        mesh->texcoords[2 * indexA + 0] = 1.0f;
        mesh->texcoords[2 * indexA + 1] = 0.0f;
        mesh->texcoords[2 * indexB + 0] = 0.0f;
        mesh->texcoords[2 * indexB + 1] = 1.0f;
        mesh->texcoords[2 * indexC + 0] = 0.0f;
        mesh->texcoords[2 * indexC + 1] = 0.0f;
    }

    UploadMesh(mesh, false);
}

static void drawQuad(Vector3 a, Vector3 b, Vector3 c, Vector3 d, Color color)
{
    DrawTriangle3D(a, b, c, color);
    DrawTriangle3D(a, c, d, color);
}

//----------------------------------------------------------------------------------
// Jugemu orbital camera
//----------------------------------------------------------------------------------
static bool moveJugemuOrbital(Camera3D *jugemu, float deltaTime)
{
    Vector3 previousPosition = jugemu->position;

    float radius           = Vector3Length(jugemu->position);
    float azimuth          = atan2f(jugemu->position.z, jugemu->position.x);
    float horizontalRadius = sqrtf(jugemu->position.x * jugemu->position.x + jugemu->position.z * jugemu->position.z);
    float elevation        = atan2f(jugemu->position.y, horizontalRadius);

    const float LONG_SPEED = 1.5f;
    const float LAT_SPEED  = 1.0f;
    const float ZOOM_SPEED = 2.0f;
    const float ROLL_SPEED = 2.0f;

    if (IsKeyDown(KEY_LEFT))
        azimuth += LONG_SPEED * deltaTime;
    if (IsKeyDown(KEY_RIGHT))
        azimuth -= LONG_SPEED * deltaTime;
    if (IsKeyDown(KEY_UP))
        elevation += LAT_SPEED * deltaTime;
    if (IsKeyDown(KEY_DOWN))
        elevation -= LAT_SPEED * deltaTime;
    if (IsKeyDown(KEY_W))
        radius -= ZOOM_SPEED * deltaTime;
    if (IsKeyDown(KEY_S))
        radius += ZOOM_SPEED * deltaTime;

    float rollDeltaRadians = 0.0f;
    if (IsKeyDown(KEY_A))
        rollDeltaRadians -= ROLL_SPEED * deltaTime;
    if (IsKeyDown(KEY_D))
        rollDeltaRadians += ROLL_SPEED * deltaTime;
    if (IsKeyPressed(KEY_SPACE))
    {
        jugemu->up       = (Vector3){0, 1, 0};
        rollDeltaRadians = 0.0f;
    }

    radius          = Clamp(radius, 0.25f, 25.0f);
    const float EPS = 0.0001f;
    elevation       = Clamp(elevation, -PI * 0.5f + EPS, PI * 0.5f - EPS);

    jugemu->position.x = radius * cosf(elevation) * cosf(azimuth);
    jugemu->position.y = radius * sinf(elevation);
    jugemu->position.z = radius * cosf(elevation) * sinf(azimuth);

    Vector3 viewDir = Vector3Normalize(Vector3Subtract((Vector3){0, 0, 0}, jugemu->position));
    Vector3 upRot   = rotatePointAboutAxis(jugemu->up, (Vector3){0, 0, 0}, viewDir, rollDeltaRadians);
    jugemu->target  = (Vector3){0, 0, 0};
    jugemu->up      = Vector3Normalize(upRot);

    return (previousPosition.x != jugemu->position.x) || (previousPosition.y != jugemu->position.y) ||
           (previousPosition.z != jugemu->position.z);
}
