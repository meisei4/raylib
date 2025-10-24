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
static const int N64_WIDTH = 640;
static const int N64_HEIGHT = 480;

static const float ANGULAR_VELOCITY = 1.25f;
static const float FOVY_PERSPECTIVE = 60.0f;

static const Vector3 MODEL_POS = {0.0f, 0.0f, 0.0f};
static const Vector3 MODEL_SCALE = {1.0f, 1.0f, 1.0f};

static const Vector3 OBSERVER_POS = {0.0f, 0.0f, 2.0f};
static const Vector3 JUGEMU_POS_ISO = {3.0f, 1.0f, 3.0f};

//----------------------------------------------------------------------------------
// Global flags/state
//----------------------------------------------------------------------------------
enum
{
    FLAG_NDC_OVERLAY = 1u << 0,         // E
    FLAG_ASPECT = 1u << 1,              // Q  (0=ISO didactic, 1=ANISO true aspect)
    FLAG_PERSPECTIVE_CORRECT = 1u << 2, // P
    FLAG_PAUSE = 1u << 3,               // F
    FLAG_COLOR_MODE = 1u << 4           // C  <-- NEW: toggle color vs texture
};

static unsigned int gFlags = 0;

#define NDC_SPACE() ((gFlags & FLAG_NDC_OVERLAY) != 0)
#define ANISOTROPIC() ((gFlags & FLAG_ASPECT) != 0)
#define PERSPECTIVE_CORRECT() ((gFlags & FLAG_PERSPECTIVE_CORRECT) != 0)
#define PAUSED() ((gFlags & FLAG_PAUSE) != 0)
#define COLOR_MODE() ((gFlags & FLAG_COLOR_MODE) != 0)
//----------------------------------------------------------------------------------
// "PerspectiveCorrect" screenshot resources CPU capture OpenGL11 + near-plane quad
//----------------------------------------------------------------------------------
static Texture2D gPerspectiveCorrectTexture = {0};

//TODO: what why? please not global shit
static Mesh gNearQuad = {0};
static Model gNearQuadModel = {0};

static bool moveJugemuOrbital(Camera3D *jugemu, float deltaTime);
static void capturePerspectiveCorrectToTexture(
    Camera3D *mainObserver,
    Model *model,
    Vector3 modelPos,
    float rotRadians,
    Vector3 modelScale,
    Texture2D currentTexture,
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
static void drawNearPlaneIntersectionalDiskMesh(
    Camera3D *observer,
    float nearClipPlane,
    const Mesh *srcMesh,
    Model *nearPlaneIntersectionModel,
    Vector3 modelPosition,
    Vector3 modelScale,
    float meshRotationRadians,
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
    float farClipPlane = 3.0f;

    Camera3D mainObserver = {0};
    mainObserver.position = OBSERVER_POS;
    mainObserver.target = (Vector3){0, 0, 0};
    mainObserver.up = (Vector3){0, 1, 0};
    mainObserver.fovy = FOVY_PERSPECTIVE;
    mainObserver.projection = CAMERA_PERSPECTIVE;

    int screenWidth = GetScreenWidth();
    int screenHeight = GetScreenHeight();
    float aspect = (float)screenWidth / (float)screenHeight;

    Camera3D jugemu = (Camera3D){0};
    jugemu.position = JUGEMU_POS_ISO;
    jugemu.target = (Vector3){0, 0, 0};
    jugemu.up = (Vector3){0, 1, 0};
    jugemu.fovy = FOVY_PERSPECTIVE;
    jugemu.projection = CAMERA_PERSPECTIVE;

    TraceLog(
        LOG_INFO,
        TextFormat("jugemu init pos: (%.3f, %.3f, %.3f)", jugemu.position.x, jugemu.position.y, jugemu.position.z));

    float idleTimer = 0.0f;
    bool movedSinceLastLog = false;

    float meshRotationRadians = 0.0f;

    Mesh cubeMesh = GenMeshCube(1.0f, 1.0f, 1.0f);
    Model mainModel = LoadModelFromMesh(cubeMesh);

    Mesh ndcMesh = (Mesh){0};
    ndcMesh.vertexCount = mainModel.meshes[0].vertexCount;
    ndcMesh.triangleCount = mainModel.meshes[0].triangleCount;
    ndcMesh.vertices = RL_CALLOC(ndcMesh.vertexCount * 3, sizeof(float));
    ndcMesh.texcoords = RL_CALLOC(ndcMesh.vertexCount * 2, sizeof(float));
    ndcMesh.indices = RL_CALLOC(ndcMesh.triangleCount * 3, sizeof(unsigned short));

    if (mainModel.meshes[0].texcoords)
    {
        for (int vertexIndex = 0; vertexIndex < ndcMesh.vertexCount; vertexIndex++)
        {
            ndcMesh.texcoords[2 * vertexIndex + 0] = mainModel.meshes[0].texcoords[2 * vertexIndex + 0];
            ndcMesh.texcoords[2 * vertexIndex + 1] = mainModel.meshes[0].texcoords[2 * vertexIndex + 1];
        }
    }
    if (mainModel.meshes[0].indices)
    {
        memcpy(ndcMesh.indices, mainModel.meshes[0].indices, sizeof(unsigned short) * 3 * ndcMesh.triangleCount);
    }

    Model ndcModel = LoadModelFromMesh(ndcMesh);

    Mesh nearPlaneIntersectionalDiskMesh = (Mesh){0};
    nearPlaneIntersectionalDiskMesh.vertexCount =
        3 * mainModel.meshes[0].triangleCount; // capacity for 3 points per tri
    nearPlaneIntersectionalDiskMesh.vertices =
        RL_CALLOC(nearPlaneIntersectionalDiskMesh.vertexCount * 3, sizeof(float));
    Model nearPlaneIntersectionalDiskModel = LoadModelFromMesh(nearPlaneIntersectionalDiskMesh);

    Image checkeredImage = GenImageChecked(16, 16, 4, 4, BLACK, WHITE);
    Texture2D checkeredTexture = LoadTextureFromImage(checkeredImage);
    UnloadImage(checkeredImage);
    mainModel.materials[0].maps[MATERIAL_MAP_ALBEDO].texture = checkeredTexture;
    ndcModel.materials[0].maps[MATERIAL_MAP_ALBEDO].texture = checkeredTexture;
    applyBarycentricPalette(&mainModel.meshes[0]);
    applyBarycentricPalette(&ndcModel.meshes[0]);
    //--------------------------------------------------------------------------------------

    // Main loop
    while (!WindowShouldClose())
    {
        // Update
        //----------------------------------------------------------------------------------
        if (IsKeyPressed(KEY_E))
            gFlags ^= FLAG_NDC_OVERLAY;
        if (IsKeyPressed(KEY_Q))
            gFlags ^= FLAG_ASPECT;
        if (IsKeyPressed(KEY_P))
            gFlags ^= FLAG_PERSPECTIVE_CORRECT;
        if (IsKeyPressed(KEY_F))
            gFlags ^= FLAG_PAUSE;
        if (IsKeyPressed(KEY_C))
            gFlags ^= FLAG_COLOR_MODE;

        screenWidth = GetScreenWidth();
        screenHeight = GetScreenHeight();
        aspect = (float)screenWidth / (float)screenHeight;

        float deltaTime = GetFrameTime();

        if (!PAUSED())
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
                TraceLog(
                    LOG_INFO,
                    TextFormat(
                        "jugemu stopped at: (%.3f, %.3f, %.3f)",
                        jugemu.position.x,
                        jugemu.position.y,
                        jugemu.position.z));
                movedSinceLastLog = false;
                idleTimer = 0.0f;
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
                &mainObserver, &mainModel, MODEL_POS, meshRotationRadians, MODEL_SCALE, checkeredTexture, COLOR_MODE());
            ClearBackground(BLACK);
        }

        BeginMode3D(jugemu);
        drawObservedAxes(&mainObserver);

        if (NDC_SPACE())
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
            //TODO: this needs work because its the whole colors vs texture and not clear on the projected near plane
            ndcModel.materials[0].maps[MATERIAL_MAP_ALBEDO].texture.id = COLOR_MODE() ? 0 : checkeredTexture.id;
            DrawModelEx(ndcModel, MODEL_POS, (Vector3){0, 1, 0}, RAD2DEG * meshRotationRadians, MODEL_SCALE, WHITE);
            rlSetLineWidth(2.0f);
            //TODO: this is not working anymore. Color mode setting the vertex array to enabled's fault???????
            ndcModel.materials[0].maps[MATERIAL_MAP_ALBEDO].texture.id = 0;
            DrawModelWiresEx(ndcModel, MODEL_POS, (Vector3){0, 1, 0}, RAD2DEG * meshRotationRadians, MODEL_SCALE, BLUE);
            rlSetPointSize(8.0f);
            DrawModelPointsEx(
                ndcModel, MODEL_POS, (Vector3){0, 1, 0}, RAD2DEG * meshRotationRadians, MODEL_SCALE, GREEN);
            drawNearPlaneIntersectionalDiskMesh(
                &mainObserver,
                nearClipPlane,
                &ndcModel.meshes[0],
                &nearPlaneIntersectionalDiskModel,
                MODEL_POS,
                MODEL_SCALE,
                meshRotationRadians,
                true);
            if (PERSPECTIVE_CORRECT())
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
                        checkeredTexture);
                }
            }
        }
        else
        {
            //TODO: this needs work because its the whole colors vs texture and not clear on the projected near plane
            mainModel.materials[0].maps[MATERIAL_MAP_ALBEDO].texture.id = COLOR_MODE() ? 0 : checkeredTexture.id;
            DrawModelEx(mainModel, MODEL_POS, (Vector3){0, 1, 0}, RAD2DEG * meshRotationRadians, MODEL_SCALE, WHITE);
            rlSetLineWidth(2.0f);
            //TODO: this is not working???? Color mode setting the vertex array to enabled's fault???????
            mainModel.materials[0].maps[MATERIAL_MAP_ALBEDO].texture.id = 0;
            DrawModelWiresEx(
                mainModel, MODEL_POS, (Vector3){0, 1, 0}, RAD2DEG * meshRotationRadians, MODEL_SCALE, BLUE);
            rlSetPointSize(8.0f);
            DrawModelPointsEx(
                mainModel, MODEL_POS, (Vector3){0, 1, 0}, RAD2DEG * meshRotationRadians, MODEL_SCALE, GREEN);

            drawFrustum(&mainObserver, ANISOTROPIC() ? aspect : 1.0f, nearClipPlane, farClipPlane);

            drawNearPlaneIntersectionalDiskMesh(
                &mainObserver,
                nearClipPlane,
                &mainModel.meshes[0],
                &nearPlaneIntersectionalDiskModel,
                MODEL_POS,
                MODEL_SCALE,
                meshRotationRadians,
                false);

            if (PERSPECTIVE_CORRECT())
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
                        checkeredTexture);
                }
            }
        }

        EndMode3D();

        // --- Key Legend ---
        DrawText(
            "E: NDC | Q: aspect | P: perspective divide | C: colors | arrowkeys/WASD: move camera | Space: reset camera",
            12,
            14,
            10,
            RAYWHITE);

        int lineY = 32;
        int lineStep = 12;

        DrawText(TextFormat("SPACE:   %s", NDC_SPACE() ? "NDC" : "WORLD"), 12, lineY, 10, RAYWHITE);
        lineY += lineStep;

        DrawText(TextFormat("ASPECT:  %s", ANISOTROPIC() ? "ANISOTROPIC" : "ISOTROPIC"), 12, lineY, 10, RAYWHITE);
        lineY += lineStep;

        DrawText(
            TextFormat("PERSPECTIVE:  %s", PERSPECTIVE_CORRECT() ? "CORRECT" : "INCORRECT"), 12, lineY, 10, RAYWHITE);
        lineY += lineStep;

        DrawText(TextFormat("COLORS:  %s", COLOR_MODE() ? "ON" : "OFF"), 12, lineY, 10, RAYWHITE);

        EndDrawing();
        //----------------------------------------------------------------------------------
    }

    // De-Initialization
    //--------------------------------------------------------------------------------------
    UnloadModel(mainModel);
    UnloadModel(ndcModel);
    UnloadModel(nearPlaneIntersectionalDiskModel);
    UnloadTexture(checkeredTexture);
    CloseWindow();
    //--------------------------------------------------------------------------------------

    return 0;
}

//----------------------------------------------------------------------------------
// Small math helpers  //TODO: many of these maybe already exist in raylib math right???
//----------------------------------------------------------------------------------
static inline Vector3 observedLineOfSight(Camera3D *observer)
{
    Vector3 lineOfSight = Vector3Subtract(observer->target, observer->position);
    return Vector3Normalize(lineOfSight);
}

static inline void observedBasis(
    Camera3D *observer, Vector3 *observedLineOfSightOut, Vector3 *observedRightOut, Vector3 *observedUpOut)
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
    //TODO: this has to exist somewhere in raylib math right???
    Vector3 axisDir = Vector3Normalize(Vector3Subtract(axisEnd, axisStart));
    Vector3 localFromAxis = Vector3Subtract(point, axisStart);
    Vector3 rotatedLocal = Vector3RotateByAxisAngle(localFromAxis, axisDir, angleRadians);
    return Vector3Add(axisStart, rotatedLocal);
}

static inline Vector3
    applyModelTranslateRotateScale(Vector3 modelCoord, Vector3 modelPosition, Vector3 modelScale, float rotationRadians)
{
    //TODO: this should be somewhere closer in proximity to where its relevant in the pipeline stages..

    Matrix M = MatrixMultiply(
        MatrixMultiply(MatrixScale(modelScale.x, modelScale.y, modelScale.z), MatrixRotateY(rotationRadians)),
        MatrixTranslate(modelPosition.x, modelPosition.y, modelPosition.z));
    return Vector3Transform(modelCoord, M);
}

static inline Vector3 applyInverseModelTranslateRotateScale(
    Vector3 worldCoord, Vector3 modelPosition, Vector3 modelScale, float rotationRadians)
{
    //TODO: this should be somewhere closer in proximity to where its relevant in the pipeline stages..

    Matrix M = MatrixMultiply(
        MatrixMultiply(MatrixScale(modelScale.x, modelScale.y, modelScale.z), MatrixRotateY(rotationRadians)),
        MatrixTranslate(modelPosition.x, modelPosition.y, modelPosition.z));
    Matrix invM = MatrixInvert(M);
    return Vector3Transform(worldCoord, invM);
}

static inline Vector3 nearPlaneIntersection(Camera3D *observer, float nearClipPlane, Vector3 worldCoord)
{
    Vector3 los = observedLineOfSight(observer);
    Vector3 rayToWorld = Vector3Subtract(worldCoord, observer->position);
    float signedDepth = Vector3DotProduct(rayToWorld, los);
    if (signedDepth <= 0.0f)
    {
        return Vector3Add(observer->position, Vector3Scale(los, nearClipPlane));
    }
    float t = nearClipPlane / signedDepth;
    return Vector3Add(observer->position, Vector3Scale(rayToWorld, t));
}

static inline Vector3 triangleNormal(Vector3 a, Vector3 b, Vector3 c)
{
    return Vector3Normalize(Vector3CrossProduct(Vector3Subtract(b, a), Vector3Subtract(c, a)));
}

static void ensureNearPlaneQuadBuilt(void)
{
    //TODO: there is likely an easier more elegant way to do this, WHY DOES THIS EVEN NEED TO EXIST...
    if (gNearQuad.vertexCount > 0 && gNearQuadModel.meshCount > 0)
        return;

    gNearQuad = (Mesh){0};
    gNearQuad.vertexCount = 4;
    gNearQuad.triangleCount = 2;
    gNearQuad.vertices = (float *)RL_CALLOC(4 * 3, sizeof(float));
    gNearQuad.texcoords = (float *)RL_CALLOC(4 * 2, sizeof(float));
    gNearQuad.indices = (unsigned short *)RL_CALLOC(6, sizeof(unsigned short));

    float uvs[8] = {0.0f, 0.0f, 1.0f, 0.0f, 1.0f, 1.0f, 0.0f, 1.0f};
    memcpy(gNearQuad.texcoords, uvs, sizeof(uvs));
    unsigned short idx[6] = {0, 2, 1, 0, 3, 2};
    memcpy(gNearQuad.indices, idx, sizeof(idx));
    gNearQuadModel = LoadModelFromMesh(gNearQuad);
}

static void updateNearPlaneQuadGeometry(Camera3D *observer, float nearClipPlane, int screenWidth, int screenHeight)
{
    //TODO: there is likely an easier more elegant way to do this, shouldnt this be gen mesh plane or something??????
    // WHY SO MUCH VECTOR MATH....

    Vector3 los, right, up;
    observedBasis(observer, &los, &right, &up);

    float aspect = (float)screenWidth / (float)screenHeight;
    float halfFovy = DEG2RAD * observer->fovy * 0.5f;
    float halfH = nearClipPlane * tanf(halfFovy);
    float halfW = ANISOTROPIC() ? (halfH * aspect) : halfH;

    Vector3 centerNear = Vector3Add(observer->position, Vector3Scale(los, nearClipPlane));

    Vector3 TL = Vector3Add(centerNear, Vector3Add(Vector3Scale(up, +halfH), Vector3Scale(right, -halfW)));
    Vector3 TR = Vector3Add(centerNear, Vector3Add(Vector3Scale(up, +halfH), Vector3Scale(right, +halfW)));
    Vector3 BR = Vector3Add(centerNear, Vector3Add(Vector3Scale(up, -halfH), Vector3Scale(right, +halfW)));
    Vector3 BL = Vector3Add(centerNear, Vector3Add(Vector3Scale(up, -halfH), Vector3Scale(right, -halfW)));

    float *v = gNearQuad.vertices;
    v[0] = TL.x;
    v[1] = TL.y;
    v[2] = TL.z;
    v[3] = TR.x;
    v[4] = TR.y;
    v[5] = TR.z;
    v[6] = BR.x;
    v[7] = BR.y;
    v[8] = BR.z;
    v[9] = BL.x;
    v[10] = BL.y;
    v[11] = BL.z;
}

static void composeAlphaFromMask(Image *rgba, const Image *mask, unsigned char threshold)
{
    if (rgba->width != mask->width || rgba->height != mask->height)
        return;

    ImageFormat(rgba, PIXELFORMAT_UNCOMPRESSED_R8G8B8A8);
    Image maskCopy = ImageCopy(*mask);
    ImageFormat(&maskCopy, PIXELFORMAT_UNCOMPRESSED_R8G8B8A8);

    unsigned char *dst = (unsigned char *)rgba->data;
    unsigned char *msk = (unsigned char *)maskCopy.data;
    int n = rgba->width * rgba->height;

    for (int i = 0; i < n; ++i)
    {
        unsigned int r = msk[4 * i + 0], g = msk[4 * i + 1], b = msk[4 * i + 2];
        unsigned int y = (r * 30u + g * 59u + b * 11u) / 100u;
        dst[4 * i + 3] = (y > threshold) ? 255 : 0;
    }

    UnloadImage(maskCopy);
}

static void capturePerspectiveCorrectToTexture(
    Camera3D *mainObserver,
    Model *model,
    Vector3 modelPos,
    float rotRadians,
    Vector3 modelScale,
    Texture2D currentTexture,
    bool usePVC)
{
    //TODO: there is likely an easier more elegant way to do this, its not intuitive to read what is happening here.

    ClearBackground(BLACK);
    BeginMode3D(*mainObserver);
    Texture2D previousTexture = model->materials[0].maps[MATERIAL_MAP_ALBEDO].texture;
    model->materials[0].maps[MATERIAL_MAP_ALBEDO].texture = usePVC ? (Texture2D){0} : currentTexture;
    DrawModelEx(*model, modelPos, (Vector3){0, 1, 0}, RAD2DEG * rotRadians, modelScale, WHITE);
    model->materials[0].maps[MATERIAL_MAP_ALBEDO].texture = previousTexture;
    EndMode3D();

    Image coloredFramebuffer = LoadImageFromScreen();
    ImageFormat(&coloredFramebuffer, PIXELFORMAT_UNCOMPRESSED_R8G8B8A8);

    ClearBackground(BLACK);
    BeginMode3D(*mainObserver);
    Texture2D cacheTexture = model->materials[0].maps[MATERIAL_MAP_ALBEDO].texture;
    Color cacheColor = model->materials[0].maps[MATERIAL_MAP_ALBEDO].color;

    model->materials[0].maps[MATERIAL_MAP_ALBEDO].texture = (Texture2D){0};
    model->materials[0].maps[MATERIAL_MAP_ALBEDO].color = WHITE;

    DrawModelEx(*model, modelPos, (Vector3){0, 1, 0}, RAD2DEG * rotRadians, modelScale, WHITE);

    model->materials[0].maps[MATERIAL_MAP_ALBEDO].texture = cacheTexture;
    model->materials[0].maps[MATERIAL_MAP_ALBEDO].color = cacheColor;
    EndMode3D();

    Image alphaMaskedFramebuffer = LoadImageFromScreen();
    composeAlphaFromMask(&coloredFramebuffer, &alphaMaskedFramebuffer, 1);

    if (gPerspectiveCorrectTexture.id != 0)
        UpdateTexture(gPerspectiveCorrectTexture, coloredFramebuffer.data);
    else
        gPerspectiveCorrectTexture = LoadTextureFromImage(coloredFramebuffer);

    UnloadImage(alphaMaskedFramebuffer);
    UnloadImage(coloredFramebuffer);
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
    //TODO: there is likely an easier more elegant way to do this, like a mesh GEN or something and modify not just drawing the lines
    Vector3 los, right, up;
    observedBasis(observer, &los, &right, &up);

    Vector3 centerNear = Vector3Add(observer->position, Vector3Scale(los, nearClipPlane));
    Vector3 centerFar = Vector3Add(observer->position, Vector3Scale(los, farClipPlane));

    float halfFovy = DEG2RAD * observer->fovy * 0.5f;
    float halfHeightNear = nearClipPlane * tanf(halfFovy);
    float halfWidthNear = halfHeightNear * aspect;
    float halfHeightFar = farClipPlane * tanf(halfFovy);
    float halfWidthFar = halfHeightFar * aspect;

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
    //TODO: there is likely an easier more elegant way to do this, like a mesh GEN and draw not just drawing the lines

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

    float halfFovy = DEG2RAD * observer->fovy * 0.5f;
    float halfHeightNear = nearClipPlane * tanf(halfFovy);
    float halfWidthNear = halfHeightNear * aspect;

    float signedDepth = Vector3DotProduct(Vector3Subtract(worldCoord, observer->position), los);
    if (signedDepth <= 0.0f)
    {
        return (Vector3){0.0f, 0.0f, 1.0f};
    }

    Vector3 intersectionCoord = nearPlaneIntersection(observer, nearClipPlane, worldCoord);
    Vector3 centerNear = Vector3Add(observer->position, Vector3Scale(los, nearClipPlane));
    Vector3 clipPlaneVector = Vector3Subtract(intersectionCoord, centerNear);

    float xNdc = Vector3DotProduct(clipPlaneVector, right) / halfWidthNear;
    float yNdc = Vector3DotProduct(clipPlaneVector, up) / halfHeightNear;

    float zNdc = ((farClipPlane + nearClipPlane) - (2.0f * farClipPlane * nearClipPlane) / signedDepth) /
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
    Vector3 los = observedLineOfSight(observer);

    float halfFovy = DEG2RAD * observer->fovy * 0.5f;
    float halfHeightNear = nearClipPlane * tanf(halfFovy);

    float halfWidthNear = ANISOTROPIC() ? (halfHeightNear * aspect) : halfHeightNear;
    float halfDepthNdc = ANISOTROPIC() ? (0.5f * (farClipPlane - nearClipPlane)) : halfHeightNear;

    Vector3 centerNear = Vector3Add(observer->position, Vector3Scale(los, nearClipPlane));
    Vector3 ndcCubeCenter = Vector3Add(centerNear, Vector3Scale(los, halfDepthNdc));
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
    float aspect = (float)screenWidth / (float)screenHeight;
    Vector3 los = observedLineOfSight(observer);
    float halfFovy = DEG2RAD * observer->fovy * 0.5f;
    float halfHeightNear = nearClipPlane * tanf(halfFovy);

    float halfWidthNear = ANISOTROPIC() ? (halfHeightNear * aspect) : halfHeightNear;
    float halfDepthNdc = ANISOTROPIC() ? (0.5f * (farClipPlane - nearClipPlane)) : halfHeightNear;

    Vector3 centerNear = Vector3Add(observer->position, Vector3Scale(los, nearClipPlane));
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

//----------------------------------------------------------------------------------
// Near-plane intersection “disk” points
//----------------------------------------------------------------------------------
static Vector3 flipYInNearClipPlane(Vector3 intersectionCoord, Camera3D *observer, float nearClipPlane)
{
    Vector3 los, right, up;
    observedBasis(observer, &los, &right, &up);
    Vector3 centerNear = Vector3Add(observer->position, Vector3Scale(los, nearClipPlane));
    Vector3 toClipPlane = Vector3Subtract(intersectionCoord, centerNear);
    float xComponent = Vector3DotProduct(toClipPlane, right);
    float yComponent = Vector3DotProduct(toClipPlane, up);
    return Vector3Add(centerNear, Vector3Add(Vector3Scale(right, xComponent), Vector3Scale(up, -yComponent)));
}
static void drawNearPlaneIntersectionalDiskMesh(
    Camera3D *observer,
    float nearClipPlane,
    const Mesh *srcMesh,
    Model *nearPlaneIntersectionModel,
    Vector3 modelPosition,
    Vector3 modelScale,
    float meshRotationRadians,
    bool reflectYAxis)
{
    Mesh *diskMesh = &nearPlaneIntersectionModel->meshes[0];
    int capacity = diskMesh->vertexCount;
    int verticesWritten = 0;

    Vector3 los = observedLineOfSight(observer);

    for (int tri = 0; tri < srcMesh->triangleCount; tri++)
    {
        int ia = srcMesh->indices ? srcMesh->indices[3 * tri + 0] : 3 * tri + 0;
        int ib = srcMesh->indices ? srcMesh->indices[3 * tri + 1] : 3 * tri + 1;
        int ic = srcMesh->indices ? srcMesh->indices[3 * tri + 2] : 3 * tri + 2;

        Vector3 A =
            (Vector3){srcMesh->vertices[3 * ia + 0], srcMesh->vertices[3 * ia + 1], srcMesh->vertices[3 * ia + 2]};
        Vector3 B =
            (Vector3){srcMesh->vertices[3 * ib + 0], srcMesh->vertices[3 * ib + 1], srcMesh->vertices[3 * ib + 2]};
        Vector3 C =
            (Vector3){srcMesh->vertices[3 * ic + 0], srcMesh->vertices[3 * ic + 1], srcMesh->vertices[3 * ic + 2]};

        A = applyModelTranslateRotateScale(A, modelPosition, modelScale, meshRotationRadians);
        B = applyModelTranslateRotateScale(B, modelPosition, modelScale, meshRotationRadians);
        C = applyModelTranslateRotateScale(C, modelPosition, modelScale, meshRotationRadians);

        Vector3 n = triangleNormal(A, B, C);
        if (Vector3DotProduct(n, los) > 0.0f)
            continue;

        Vector3 hits[3] = {
            nearPlaneIntersection(observer, nearClipPlane, A),
            nearPlaneIntersection(observer, nearClipPlane, B),
            nearPlaneIntersection(observer, nearClipPlane, C)};

        for (int v = 0; v < 3 && verticesWritten < capacity; ++v)
        {
            Vector3 worldV = (Vector3[]){A, B, C}[v];
            Vector3 hit = reflectYAxis ? flipYInNearClipPlane(hits[v], observer, nearClipPlane) : hits[v];

            rlSetLineWidth(1.0f);
            DrawLine3D(worldV, hit, (Color){255, 0, 0, 80});

            diskMesh->vertices[3 * verticesWritten + 0] = hit.x;
            diskMesh->vertices[3 * verticesWritten + 1] = hit.y;
            diskMesh->vertices[3 * verticesWritten + 2] = hit.z;
            verticesWritten++;
        }
    }

    diskMesh->vertexCount = verticesWritten;

    rlSetPointSize(6.0f);
    DrawModelPoints(*nearPlaneIntersectionModel, (Vector3){0, 0, 0}, 1.0f, GREEN);
}

static inline Vector3 remapNearPlaneByAspect(Vector3 hit, Vector3 centerNear, Vector3 right, Vector3 up, float xScale)
{
    //TODO: this is gross, and naming of functions is not clear and this is just bad and not clear at all in the chornology of the pipeline didactic
    Vector3 d = Vector3Subtract(hit, centerNear);
    float xr = Vector3DotProduct(d, right);
    float yu = Vector3DotProduct(d, up);
    return Vector3Add(centerNear, Vector3Add(Vector3Scale(right, xr * xScale), Vector3Scale(up, yu)));
}

//TODO: this LITERALLY does not work. color UVs and TEXCOORDS are NOT THE SAME (this is a good demonstration of why that is so i think????
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
        applyBarycentricPalette(
            mesh); //TODO: remove this and just be more clear about toggling behavior idempotency or something

    float aspect = (float)screenWidth / (float)screenHeight;
    Vector3 los, right, up;
    observedBasis(observer, &los, &right, &up);

    Vector3 centerNear = Vector3Add(observer->position, Vector3Scale(los, nearClipPlane));
    float xScale = ANISOTROPIC() ? 1.0f : (1.0f / aspect);

    rlDisableTexture();
    rlBegin(RL_TRIANGLES);

    for (int tri = 0; tri < mesh->triangleCount; tri++)
    {
        int ia = mesh->indices ? mesh->indices[3 * tri + 0] : 3 * tri + 0;
        int ib = mesh->indices ? mesh->indices[3 * tri + 1] : 3 * tri + 1;
        int ic = mesh->indices ? mesh->indices[3 * tri + 2] : 3 * tri + 2;

        Vector3 A = (Vector3){mesh->vertices[3 * ia + 0], mesh->vertices[3 * ia + 1], mesh->vertices[3 * ia + 2]};
        Vector3 B = (Vector3){mesh->vertices[3 * ib + 0], mesh->vertices[3 * ib + 1], mesh->vertices[3 * ib + 2]};
        Vector3 C = (Vector3){mesh->vertices[3 * ic + 0], mesh->vertices[3 * ic + 1], mesh->vertices[3 * ic + 2]};

        A = applyModelTranslateRotateScale(A, modelPosition, modelScale, meshRotationRadians);
        B = applyModelTranslateRotateScale(B, modelPosition, modelScale, meshRotationRadians);
        C = applyModelTranslateRotateScale(C, modelPosition, modelScale, meshRotationRadians);

        Vector3 hitA = nearPlaneIntersection(observer, nearClipPlane, A);
        Vector3 hitB = nearPlaneIntersection(observer, nearClipPlane, B);
        Vector3 hitC = nearPlaneIntersection(observer, nearClipPlane, C);

        Vector3 pA = remapNearPlaneByAspect(hitA, centerNear, right, up, xScale);
        Vector3 pB = remapNearPlaneByAspect(hitB, centerNear, right, up, xScale);
        Vector3 pC = remapNearPlaneByAspect(hitC, centerNear, right, up, xScale);

        const unsigned char *cA = &mesh->colors[4 * ia];
        const unsigned char *cB = &mesh->colors[4 * ib];
        const unsigned char *cC = &mesh->colors[4 * ic];

        rlColor4ub(cA[0], cA[1], cA[2], cA[3]);
        rlVertex3f(pA.x, pA.y, pA.z);
        rlColor4ub(cB[0], cB[1], cB[2], cB[3]);
        rlVertex3f(pB.x, pB.y, pB.z);
        rlColor4ub(cC[0], cC[1], cC[2], cC[3]);
        rlVertex3f(pC.x, pC.y, pC.z);
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

    Vector3 centerNear = Vector3Add(observer->position, Vector3Scale(los, nearClipPlane));
    float xScale = ANISOTROPIC() ? 1.0f : (1.0f / aspect);

    rlColor4f(1, 1, 1, 1);
    rlEnableTexture(texture.id);
    rlSetTexture(texture.id);
    rlBegin(RL_TRIANGLES);

    for (int triangleIndex = 0; triangleIndex < mesh->triangleCount; triangleIndex++)
    {
        int indexA = mesh->indices ? mesh->indices[3 * triangleIndex + 0] : 3 * triangleIndex + 0;
        int indexB = mesh->indices ? mesh->indices[3 * triangleIndex + 1] : 3 * triangleIndex + 1;
        int indexC = mesh->indices ? mesh->indices[3 * triangleIndex + 2] : 3 * triangleIndex + 2;

        Vector3 vertexA =
            (Vector3){mesh->vertices[3 * indexA + 0], mesh->vertices[3 * indexA + 1], mesh->vertices[3 * indexA + 2]};
        Vector3 vertexB =
            (Vector3){mesh->vertices[3 * indexB + 0], mesh->vertices[3 * indexB + 1], mesh->vertices[3 * indexB + 2]};
        Vector3 vertexC =
            (Vector3){mesh->vertices[3 * indexC + 0], mesh->vertices[3 * indexC + 1], mesh->vertices[3 * indexC + 2]};

        vertexA = applyModelTranslateRotateScale(vertexA, modelPosition, modelScale, meshRotationRadians);
        vertexB = applyModelTranslateRotateScale(vertexB, modelPosition, modelScale, meshRotationRadians);
        vertexC = applyModelTranslateRotateScale(vertexC, modelPosition, modelScale, meshRotationRadians);

        Vector3 hitA = nearPlaneIntersection(observer, nearClipPlane, vertexA);
        Vector3 hitB = nearPlaneIntersection(observer, nearClipPlane, vertexB);
        Vector3 hitC = nearPlaneIntersection(observer, nearClipPlane, vertexC);

        Vector3 planeApos = remapNearPlaneByAspect(hitA, centerNear, right, up, xScale);
        Vector3 planeBpos = remapNearPlaneByAspect(hitB, centerNear, right, up, xScale);
        Vector3 planeCpos = remapNearPlaneByAspect(hitC, centerNear, right, up, xScale);

        float uA = mesh->texcoords ? mesh->texcoords[2 * indexA + 0] : 0.0f;
        float vA = mesh->texcoords ? mesh->texcoords[2 * indexA + 1] : 0.0f;
        float uB = mesh->texcoords ? mesh->texcoords[2 * indexB + 0] : 0.0f;
        float vB = mesh->texcoords ? mesh->texcoords[2 * indexB + 1] : 0.0f;
        float uC = mesh->texcoords ? mesh->texcoords[2 * indexC + 0] : 0.0f;
        float vC = mesh->texcoords ? mesh->texcoords[2 * indexC + 1] : 0.0f;

        rlTexCoord2f(uA, vA);
        rlVertex3f(planeApos.x, planeApos.y, planeApos.z);
        rlTexCoord2f(uB, vB);
        rlVertex3f(planeBpos.x, planeBpos.y, planeBpos.z);
        rlTexCoord2f(uC, vC);
        rlVertex3f(planeCpos.x, planeCpos.y, planeCpos.z);
    }

    rlEnd();
    rlDisableTexture();
}

//----------------------------------------------------------------------------------
// Palette coloring by triangle; initializes colors and texcoords if needed
//----------------------------------------------------------------------------------
static void applyBarycentricPalette(Mesh *mesh)
{
    if (!mesh || mesh->vertexCount <= 0 || mesh->triangleCount <= 0)
        return;

    if (mesh->colors == NULL)
        mesh->colors = (unsigned char *)RL_CALLOC(mesh->vertexCount * 4, sizeof(unsigned char));
    //TODO: fix all shorthand variable naming throughout the whole code, this is not good, and not consistent please
    for (int tri = 0; tri < mesh->triangleCount; tri++)
    {
        int ia = mesh->indices ? mesh->indices[3 * tri + 0] : 3 * tri + 0;
        int ib = mesh->indices ? mesh->indices[3 * tri + 1] : 3 * tri + 1;
        int ic = mesh->indices ? mesh->indices[3 * tri + 2] : 3 * tri + 2;

        //TODO: can we more cleanly achieve this? like i just want to set a color like = RED = GREEN = BLUE, i dont want ot fuck with pointer offsets and arrays...
        // A -> Red
        mesh->colors[4 * ia + 0] = 255;
        mesh->colors[4 * ia + 1] = 0;
        mesh->colors[4 * ia + 2] = 0;
        mesh->colors[4 * ia + 3] = 255;

        // B -> Green
        mesh->colors[4 * ib + 0] = 0;
        mesh->colors[4 * ib + 1] = 255;
        mesh->colors[4 * ib + 2] = 0;
        mesh->colors[4 * ib + 3] = 255;

        // C -> Blue
        mesh->colors[4 * ic + 0] = 0;
        mesh->colors[4 * ic + 1] = 0;
        mesh->colors[4 * ic + 2] = 255;
        mesh->colors[4 * ic + 3] = 255;
    }
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
    const float LAT_SPEED = 1.0f;
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
        jugemu->up = (Vector3){0, 1, 0};
        rollDeltaRadians = 0.0f;
    }

    radius = Clamp(radius, 0.25f, 25.0f);
    const float EPS = 0.0001f;
    elevation = Clamp(elevation, -PI * 0.5f + EPS, PI * 0.5f - EPS);

    jugemu->position.x = radius * cosf(elevation) * cosf(azimuth);
    jugemu->position.y = radius * sinf(elevation);
    jugemu->position.z = radius * cosf(elevation) * sinf(azimuth);

    Vector3 viewDir = Vector3Normalize(Vector3Subtract((Vector3){0, 0, 0}, jugemu->position));
    Vector3 upRot = rotatePointAboutAxis(jugemu->up, (Vector3){0, 0, 0}, viewDir, rollDeltaRadians);
    jugemu->target = (Vector3){0, 0, 0};
    jugemu->up = Vector3Normalize(upRot);

    return (previousPosition.x != jugemu->position.x) || (previousPosition.y != jugemu->position.y) ||
           (previousPosition.z != jugemu->position.z);
}