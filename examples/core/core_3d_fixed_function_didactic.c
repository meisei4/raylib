#include <float.h>

#include "raylib.h"
#include "raymath.h"
#include "rlgl.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

//TODO:
// 0. REWRITE ALL THE POINTER ARITHMETIC WITH ACTUAL TYPES AND .x, .y., .z, .r .g. b etc etc
// 1. add proper clipping to the target meshes to show how it works with the spaces
//   - move mesh out of clip planes or allow moving the main camera's target away from the meshes
// 2. improve didactic annotations (ideally with spatial labeling rather than simple flat screen overlay)
// 3. improve code didactic, code should read in order of fixed function staging
// 4. add scripted toggling/navigation of ordered fixed function staging visualization
// 5. add some sort of ghosting effect between fixed function stages, to emphasize previous stage perhaps)
// 6. OPTIONALLY improve toggling and space navigation

static const int FONT_SIZE = 20;
static const int WIDTH = 800;
static const int HEIGHT = 450;
static const float ANGULAR_VELOCITY = 1.25f;
static const float FOVY_PERSPECTIVE = 60.0f;
static const Vector3 MODEL_POS = {0.0f, 0.0f, 0.0f};
static const Vector3 MODEL_SCALE = {1.0f, 1.0f, 1.0f};
static const Vector3 MAIN_POS = {0.0f, 0.0f, 2.0f};
static const Vector3 JUGEMU_POS_ISO = {3.0f, 1.0f, 3.0f};
static const float BLEND_SCALAR = 5.0f;

enum
{
    FLAG_NDC = 1u << 0,
    FLAG_REFLECT_Y = 1u << 1,
    FLAG_ASPECT = 1u << 2,
    FLAG_PERSPECTIVE_CORRECT = 1u << 3,
    FLAG_PAUSE = 1u << 4,
    FLAG_COLOR_MODE = 1u << 5,
    FLAG_TEXTURE_MODE = 1u << 6
};
static unsigned int gFlags = FLAG_ASPECT | FLAG_COLOR_MODE;

#define NDC_SPACE() ((gFlags & FLAG_NDC) != 0)
#define REFLECT_Y() ((gFlags & FLAG_REFLECT_Y) != 0)
#define ASPECT_CORRECT() ((gFlags & FLAG_ASPECT) != 0)
#define PERSPECTIVE_CORRECT() ((gFlags & FLAG_PERSPECTIVE_CORRECT) != 0)
#define PAUSED() ((gFlags & FLAG_PAUSE) != 0)
#define COLOR_MODE() ((gFlags & FLAG_COLOR_MODE) != 0)
#define TEXTURE_MODE() ((gFlags & FLAG_TEXTURE_MODE) != 0)
#define TOGGLE(K, F)                                                                                                                                           \
    do                                                                                                                                                         \
    {                                                                                                                                                          \
        if (IsKeyPressed(K)) gFlags ^= (F);                                                                                                                    \
    } while (0)

#define BAHAMA_BLUE CLITERAL(Color){0, 102, 153, 255}
#define SUNFLOWER CLITERAL(Color){255, 204, 153, 255}
#define PALE_CANARY CLITERAL(Color){255, 255, 153, 255}
#define ANAKIWA CLITERAL(Color){153, 204, 255, 255}
#define MARINER CLITERAL(Color){51, 102, 204, 255}
#define NEON_CARROT CLITERAL(Color){255, 153, 51, 255}
#define EGGPLANT CLITERAL(Color){102, 68, 102, 255}
#define HOPBUSH CLITERAL(Color){204, 102, 153, 255}
#define LILAC CLITERAL(Color){204, 153, 204, 255}
#define RED_DAMASK CLITERAL(Color){221, 102, 68, 255}
#define CHESTNUT_ROSE CLITERAL(Color){204, 102, 102, 255}

static void move(Camera3D *jugemu, float deltaTime);
static void drawBasisVector(const Camera3D *main);
static void worldToNDCSpace(const Camera3D *main, float aspect, float near, float far, const Model *world, const Model *ndc, float rotation);
static void drawModelFilled(const Model *model, Texture2D texture, float rotation);
static void drawModelWiresAndPoints(const Model *model, float rotation);
static void drawNearPlanePoints(const Camera3D *main, float aspect, float near, const Model *nearPlanePointsModel, const Mesh *mesh, float rotation);
static void updateSpatialFrame(const Camera3D *main, float aspect, float near, float far, const Mesh *spatialFrame);
static void drawSpatialFrame(Mesh spatialFrame);
static void drawNearPlane(
    const Camera3D *main,
    float aspect,
    float near,
    Model spatialFrameModel,
    const Mesh *mesh,
    Texture2D meshTexture,
    Texture2D perspectiveCorrectTexture,
    float rotation);
static void perspectiveIncorrectCapture(const Camera3D *main, float aspect, float near, const Mesh *mesh, Texture2D meshTexture, float rotation);
static void perspectiveCorrectCapture(const Camera3D *main, const Model *model, Texture2D meshTexture, Texture2D *perspectiveCorrectTexture, float rotation);
static void alphaMaskPunchOut(Image *rgba, const Image *mask, unsigned char threshold);
static void fillVertexColors(Mesh *mesh);
static float spaceBlendFactor(float deltaTime);
static float aspectBlendFactor(float deltaTime);
static float reflectBlendFactor(float deltaTime);

int main(void)
{
    InitWindow(WIDTH, HEIGHT, "fixed function didactic");
    SetTargetFPS(60);

    Camera3D main = {0};
    main.position = MAIN_POS;
    main.target = (Vector3){0, 0, 0};
    main.up = (Vector3){0, 1, 0};
    main.fovy = FOVY_PERSPECTIVE;
    main.projection = CAMERA_PERSPECTIVE;
    float aspect = (float)GetScreenWidth() / (float)GetScreenHeight();
    Camera3D jugemu = (Camera3D){0};
    jugemu.position = JUGEMU_POS_ISO;
    jugemu.target = (Vector3){0, 0, 0};
    jugemu.up = (Vector3){0, 1, 0};
    jugemu.fovy = FOVY_PERSPECTIVE;
    jugemu.projection = CAMERA_PERSPECTIVE;
    float meshRotation = 0.0f;
    // Model worldModel = LoadModelFromMesh(GenMeshCube(1.0f, 1.0f, 1.0f));
    // Image textureImage = GenImageChecked(4, 4, 1, 1, BLACK, WHITE);

    // Model worldModel = LoadModelFromMesh(GenMeshSphere(0.5, 8, 8));
    // Image textureImage = GenImageChecked(16, 16, 1, 1, BLACK, WHITE);

    // Model worldModel = LoadModelFromMesh(GenMeshKnot(1.0f, 1.0f, 8, 64));
    Model worldModel = LoadModelFromMesh(GenMeshKnot(1.0f, 1.0f, 16, 128));
    Image textureImage = GenImageChecked(32, 32, 1, 1, BLACK, WHITE);

    Texture2D meshTexture = LoadTextureFromImage(textureImage);
    UnloadImage(textureImage);

    Mesh ndcMesh = (Mesh){0};
    ndcMesh.vertexCount = worldModel.meshes[0].vertexCount;
    ndcMesh.triangleCount = worldModel.meshes[0].triangleCount;
    ndcMesh.vertices = RL_CALLOC(ndcMesh.vertexCount * 3, sizeof(float));
    ndcMesh.texcoords = RL_CALLOC(ndcMesh.vertexCount * 2, sizeof(float));
    ndcMesh.indices = RL_CALLOC(ndcMesh.triangleCount * 3, sizeof(unsigned short));
    ndcMesh.colors = (unsigned char *)RL_CALLOC(ndcMesh.vertexCount, sizeof(Color));
    if (worldModel.meshes[0].texcoords) memcpy(ndcMesh.texcoords, worldModel.meshes[0].texcoords, sizeof(float) * 2 * ndcMesh.vertexCount);
    if (worldModel.meshes[0].indices)
        memcpy(ndcMesh.indices, worldModel.meshes[0].indices, sizeof(unsigned short) * 3 * ndcMesh.triangleCount);
    else
        ndcMesh.indices = NULL;
    fillVertexColors(&worldModel.meshes[0]);
    memcpy(ndcMesh.colors, worldModel.meshes[0].colors, ndcMesh.vertexCount * sizeof(Color));

    Model ndcModel = LoadModelFromMesh(ndcMesh);

    Mesh nearPlanePoints = (Mesh){0};
    nearPlanePoints.vertexCount = 3 * worldModel.meshes[0].triangleCount;
    nearPlanePoints.vertices = RL_CALLOC(nearPlanePoints.vertexCount * 3, sizeof(float));
    Model nearPlanePointsModel = LoadModelFromMesh(nearPlanePoints);

    worldModel.materials[0].maps[MATERIAL_MAP_ALBEDO].texture = meshTexture;
    ndcModel.materials[0].maps[MATERIAL_MAP_ALBEDO].texture = meshTexture;

    Texture2D perspectiveCorrectTexture = {0}; //TODO: still wacky to do this and pass by reference everywhere i dont like that
    Mesh spatialFrame = GenMeshCube(1.0f, 1.0f, 1.0f);
    spatialFrame.colors = (unsigned char *)RL_CALLOC(spatialFrame.vertexCount, sizeof(Color));
    for (int i = 0; i < spatialFrame.vertexCount; i++) ((Color *)spatialFrame.colors)[i] = (Color){255, 255, 255, 0};
    for (int i = 0; i < 4; i++) ((Color *)spatialFrame.colors)[i].a = 255;
    Model spatialFrameModel = LoadModelFromMesh(spatialFrame);
    while (!WindowShouldClose())
    {
        const float far = 3.0f;
        const float near = 1.0f;
        aspect = (float)GetScreenWidth() / (float)GetScreenHeight();
        TOGGLE(KEY_N, FLAG_NDC);
        if (NDC_SPACE()) TOGGLE(KEY_F, FLAG_REFLECT_Y);
        TOGGLE(KEY_Q, FLAG_ASPECT);
        TOGGLE(KEY_P, FLAG_PERSPECTIVE_CORRECT);
        TOGGLE(KEY_SPACE, FLAG_PAUSE);
        TOGGLE(KEY_C, FLAG_COLOR_MODE);
        TOGGLE(KEY_T, FLAG_TEXTURE_MODE);

        move(&jugemu, GetFrameTime());
        if (!PAUSED()) meshRotation -= ANGULAR_VELOCITY * GetFrameTime();
        worldToNDCSpace(&main, aspect, near, far, &worldModel, &ndcModel, meshRotation);
        float blend = spaceBlendFactor(GetFrameTime());
        aspectBlendFactor(GetFrameTime());
        reflectBlendFactor(GetFrameTime());
        for (int i = 0; i < ndcModel.meshes[0].vertexCount; ++i)
        {
            float worldX = worldModel.meshes[0].vertices[3 * i + 0], worldY = worldModel.meshes[0].vertices[3 * i + 1],
                  worldZ = worldModel.meshes[0].vertices[3 * i + 2];
            float ndcX = ndcModel.meshes[0].vertices[3 * i + 0], ndcY = ndcModel.meshes[0].vertices[3 * i + 1], ndcZ = ndcModel.meshes[0].vertices[3 * i + 2];
            ndcModel.meshes[0].vertices[3 * i + 0] = Lerp(worldX, ndcX, blend);
            ndcModel.meshes[0].vertices[3 * i + 1] = Lerp(worldY, ndcY, blend);
            ndcModel.meshes[0].vertices[3 * i + 2] = Lerp(worldZ, ndcZ, blend);
        }
        Model *displayModel = &ndcModel;
        Mesh *displayMesh = &ndcModel.meshes[0];

        BeginDrawing();
        if (PERSPECTIVE_CORRECT() && TEXTURE_MODE()) perspectiveCorrectCapture(&main, displayModel, meshTexture, &perspectiveCorrectTexture, meshRotation);
        ClearBackground(BLACK);

        BeginMode3D(jugemu);
        drawBasisVector(&main);
        updateSpatialFrame(&main, aspect, near, far, &spatialFrame);
        drawSpatialFrame(spatialFrame);
        drawModelFilled(displayModel, meshTexture, meshRotation);
        drawModelWiresAndPoints(displayModel, meshRotation);
        drawNearPlanePoints(&main, aspect, near, &nearPlanePointsModel, displayMesh, meshRotation);
        drawNearPlane(&main, aspect, near, spatialFrameModel, displayMesh, meshTexture, perspectiveCorrectTexture, meshRotation);
        EndMode3D();

        DrawText("ARROWS: MOVE | SPACEBAR: PAUSE", 12, 12, FONT_SIZE, NEON_CARROT);
        DrawText("W A : ZOOM", 12, 38, FONT_SIZE, NEON_CARROT);
        DrawText("TEXTURE [ T ]:", 570, 12, FONT_SIZE, SUNFLOWER);
        DrawText(TEXTURE_MODE() ? "ON" : "OFF", 740, 12, FONT_SIZE, TEXTURE_MODE() ? ANAKIWA : CHESTNUT_ROSE);
        DrawText("COLORS [ C ]:", 570, 38, FONT_SIZE, SUNFLOWER);
        DrawText(COLOR_MODE() ? "ON" : "OFF", 740, 38, FONT_SIZE, COLOR_MODE() ? ANAKIWA : CHESTNUT_ROSE);
        DrawText("ASPECT [ Q ]:", 12, 392, FONT_SIZE, SUNFLOWER);
        DrawText(ASPECT_CORRECT() ? "CORRECT" : "INCORRECT", 230, 392, FONT_SIZE, ASPECT_CORRECT() ? ANAKIWA : CHESTNUT_ROSE);
        DrawText("PERSPECTIVE [ P ]:", 12, 418, FONT_SIZE, SUNFLOWER);
        DrawText(PERSPECTIVE_CORRECT() ? "CORRECT" : "INCORRECT", 230, 418, FONT_SIZE, PERSPECTIVE_CORRECT() ? ANAKIWA : CHESTNUT_ROSE);
        DrawText("SPACE [ N ]:", 530, 392, FONT_SIZE, SUNFLOWER);
        DrawText(NDC_SPACE() ? "NDC" : "WORLD", 665, 392, FONT_SIZE, NDC_SPACE() ? BAHAMA_BLUE : ANAKIWA);
        if (NDC_SPACE())
        {
            DrawText("REFLECT [ F ]:", 530, 418, FONT_SIZE, SUNFLOWER);
            DrawText(REFLECT_Y() ? "Y_DOWN" : "Y_UP", 695, 418, FONT_SIZE, REFLECT_Y() ? ANAKIWA : CHESTNUT_ROSE);
        }
        EndDrawing();
    }

    UnloadModel(worldModel);
    UnloadModel(ndcModel);
    UnloadModel(nearPlanePointsModel);
    UnloadModel(spatialFrameModel);
    UnloadTexture(meshTexture);
    CloseWindow();
    return 0;
}

static void move(Camera3D *jugemu, const float deltaTime)
{
    float radius = Vector3Length(jugemu->position);
    float azimuth = atan2f(jugemu->position.z, jugemu->position.x);
    const float horizontalRadius = sqrtf(jugemu->position.x * jugemu->position.x + jugemu->position.z * jugemu->position.z);
    float elevation = atan2f(jugemu->position.y, horizontalRadius);
    const float LONG_SPEED = 1.5f;
    const float LAT_SPEED = 1.0f;
    const float ZOOM_SPEED = 2.0f;
    if (IsKeyDown(KEY_LEFT)) azimuth += LONG_SPEED * deltaTime;
    if (IsKeyDown(KEY_RIGHT)) azimuth -= LONG_SPEED * deltaTime;
    if (IsKeyDown(KEY_UP)) elevation += LAT_SPEED * deltaTime;
    if (IsKeyDown(KEY_DOWN)) elevation -= LAT_SPEED * deltaTime;
    if (IsKeyDown(KEY_W)) radius -= ZOOM_SPEED * deltaTime;
    if (IsKeyDown(KEY_S)) radius += ZOOM_SPEED * deltaTime;
    elevation = Clamp(elevation, -M_PI_2 + 0.1f, M_PI_2 - 0.1f);
    jugemu->position.x = Clamp(radius, 0.25f, 10.0f) * cosf(elevation) * cosf(azimuth);
    jugemu->position.y = Clamp(radius, 0.25f, 10.0f) * sinf(elevation);
    jugemu->position.z = Clamp(radius, 0.25f, 10.0f) * cosf(elevation) * sinf(azimuth);
}

static void basisVector(const Camera3D *main, Vector3 *depthOut, Vector3 *rightOut, Vector3 *upOut)
{
    const Vector3 depth = Vector3Normalize(Vector3Subtract(main->target, main->position));
    const Vector3 right = Vector3Normalize(Vector3CrossProduct(depth, main->up));
    const Vector3 up = Vector3Normalize(Vector3CrossProduct(right, depth));
    *depthOut = depth;
    *rightOut = right;
    *upOut = up;
}

static void drawBasisVector(const Camera3D *main)
{
    Vector3 depth, right, up;
    basisVector(main, &depth, &right, &up);
    DrawLine3D(main->position, Vector3Add(main->position, right), NEON_CARROT);
    DrawLine3D(main->position, Vector3Add(main->position, up), LILAC);
    DrawLine3D(main->position, Vector3Add(main->position, depth), MARINER);
}

static Vector3 translateRotateScale(const int inverse, const Vector3 coordinate, const Vector3 position, const Vector3 scale, const float rotation)
{
    const Matrix matrix =
        MatrixMultiply(MatrixMultiply(MatrixScale(scale.x, scale.y, scale.z), MatrixRotateY(rotation)), MatrixTranslate(position.x, position.y, position.z));
    const Matrix result = inverse ? MatrixInvert(matrix) : matrix;
    return Vector3Transform(coordinate, result);
}

static Vector3 nearPlaneIntersection(const Camera3D *main, const float near, const Vector3 worldCoord)
{
    const Vector3 viewDir = Vector3Normalize(Vector3Subtract(main->target, main->position));
    const Vector3 mainCameraToPoint = Vector3Subtract(worldCoord, main->position);
    const float depthAlongView = Vector3DotProduct(mainCameraToPoint, viewDir);
    if (depthAlongView <= 0.0f) return Vector3Add(main->position, Vector3Scale(viewDir, near));
    const float scaleToNear = near / depthAlongView;
    return Vector3Add(main->position, Vector3Scale(mainCameraToPoint, scaleToNear));
}

static void
    worldToNDCSpace(const Camera3D *main, const float aspect, const float near, const float far, const Model *world, const Model *ndc, const float rotation)
{
    Vector3 depth, right, up;
    basisVector(main, &depth, &right, &up);
    const float halfHNear = near * tanf(DEG2RAD * main->fovy * 0.5f);
    const float halfWNear = Lerp(halfHNear, halfHNear * aspect, aspectBlendFactor(0.0f));
    const float halfDepthNDC = Lerp(halfHNear, 0.5f * (far - near), aspectBlendFactor(0.0f));
    const Vector3 centerNearPlane = Vector3Add(main->position, Vector3Scale(depth, near));
    const Vector3 centerNDCCube = Vector3Add(centerNearPlane, Vector3Scale(depth, halfDepthNDC));
    const Mesh *worldMesh = &world->meshes[0];
    const Mesh *ndcMesh = &ndc->meshes[0];

    for (int i = 0; i < worldMesh->vertexCount; i++)
    {
        const Vector3 objectVertex = {worldMesh->vertices[3 * i + 0], worldMesh->vertices[3 * i + 1], worldMesh->vertices[3 * i + 2]};
        const Vector3 worldVertex = translateRotateScale(0, objectVertex, MODEL_POS, MODEL_SCALE, rotation);
        const float signedDepth = Vector3DotProduct(Vector3Subtract(worldVertex, main->position), depth);
        const Vector3 intersectionCoord = nearPlaneIntersection(main, near, worldVertex);
        const Vector3 clipPlaneVector = Vector3Subtract(intersectionCoord, centerNearPlane);
        const float xNDC = Vector3DotProduct(clipPlaneVector, right) / halfWNear;
        const float yNDC = Vector3DotProduct(clipPlaneVector, up) / halfHNear;
        const float zNDC = ((far + near) - (2.0f * far * near) / signedDepth) / (far - near);
        const Vector3 ndcCoord = (Vector3){xNDC, yNDC, zNDC};

        const Vector3 scaledRight = Vector3Scale(right, ndcCoord.x * halfWNear);
        const Vector3 scaledUp = Vector3Scale(up, ndcCoord.y * halfHNear);
        const Vector3 scaledDepth = Vector3Scale(depth, ndcCoord.z * halfDepthNDC);
        const Vector3 offset = Vector3Add(Vector3Add(scaledRight, scaledUp), scaledDepth);
        const Vector3 scaledNDCCoord = Vector3Add(centerNDCCube, offset);

        const Vector3 mappedObjectCoord = translateRotateScale(1, scaledNDCCoord, MODEL_POS, MODEL_SCALE, rotation);
        ndcMesh->vertices[3 * i + 0] = mappedObjectCoord.x;
        ndcMesh->vertices[3 * i + 1] = mappedObjectCoord.y;
        ndcMesh->vertices[3 * i + 2] = mappedObjectCoord.z;
    }
}

static void drawModelFilled(const Model *model, const Texture2D texture, const float rotation)
{
    if (!(COLOR_MODE() || TEXTURE_MODE())) return;
    unsigned char *colorsBackup = model->meshes[0].colors;
    if (TEXTURE_MODE() && !COLOR_MODE()) model->meshes[0].colors = NULL;
    model->materials[0].maps[MATERIAL_MAP_ALBEDO].texture.id = TEXTURE_MODE() ? texture.id : 0;
    DrawModelEx(*model, MODEL_POS, (Vector3){0, 1, 0}, RAD2DEG * rotation, MODEL_SCALE, WHITE);
    model->materials[0].maps[MATERIAL_MAP_ALBEDO].texture.id = 0;
    model->meshes[0].colors = colorsBackup;
}

static void drawModelWiresAndPoints(const Model *model, const float rotation)
{
    unsigned char *cacheColors = model->meshes[0].colors;
    const unsigned int cacheID = model->materials[0].maps[MATERIAL_MAP_ALBEDO].texture.id;
    model->meshes[0].colors = NULL;
    model->materials[0].maps[MATERIAL_MAP_ALBEDO].texture.id = 0;
    DrawModelWiresEx(*model, MODEL_POS, (Vector3){0, 1, 0}, RAD2DEG * rotation, MODEL_SCALE, MARINER);
    rlSetPointSize(4.0f);
    DrawModelPointsEx(*model, MODEL_POS, (Vector3){0, 1, 0}, RAD2DEG * rotation, MODEL_SCALE, LILAC);
    model->meshes[0].colors = cacheColors;
    model->materials[0].maps[MATERIAL_MAP_ALBEDO].texture.id = cacheID;
}

static void updateSpatialFrame(const Camera3D *main, const float aspect, const float near, const float far, const Mesh *spatialFrame)
{
    Vector3 depth, right, up;
    basisVector(main, &depth, &right, &up);
    const float halfHNear = near * tanf(DEG2RAD * main->fovy * 0.5f);
    const float halfHFar = far * tanf(DEG2RAD * main->fovy * 0.5f);

    const float halfWNear = Lerp(halfHNear, halfHNear * aspect, aspectBlendFactor(0.0f));
    const float halfWFar = Lerp(halfHFar, halfHFar * aspect, aspectBlendFactor(0.0f));

    const Vector3 centerNear = Vector3Add(main->position, Vector3Scale(depth, near));
    const float halfDepthWorld = 0.5f * (far - near);
    const float halfDepthNdc = Lerp(halfHNear, 0.5f * (far - near), aspectBlendFactor(0.0f));
    const float halfDepth = Lerp(halfDepthWorld, halfDepthNdc, spaceBlendFactor(0.0f));

    const float depthSpan = 2.0f * halfDepth;
    const float farHalfW = Lerp(halfWFar, halfWNear, spaceBlendFactor(0.0f));
    const float farHalfH = Lerp(halfHFar, halfHNear, spaceBlendFactor(0.0f));

    for (int i = 0; i < spatialFrame->vertexCount; ++i)
    {
        const Vector3 vertex = (Vector3){spatialFrame->vertices[3 * i + 0], spatialFrame->vertices[3 * i + 1], spatialFrame->vertices[3 * i + 2]};
        const Vector3 offset = Vector3Subtract(vertex, centerNear);
        const float xSign = (Vector3DotProduct(offset, right) >= 0.0f) ? 1.0f : -1.0f;
        const float ySign = (Vector3DotProduct(offset, up) >= 0.0f) ? 1.0f : -1.0f;
        const float farMask = (Vector3DotProduct(offset, depth) > halfDepth) ? 1.0f : 0.0f;
        const float finalHalfWidth = halfWNear + farMask * (farHalfW - halfWNear);
        const float finalHalfHeight = halfHNear + farMask * (farHalfH - halfHNear);
        const Vector3 faceCenter = Vector3Add(centerNear, Vector3Scale(depth, farMask * depthSpan));
        const Vector3 result = Vector3Add(faceCenter, Vector3Add(Vector3Scale(right, xSign * finalHalfWidth), Vector3Scale(up, ySign * finalHalfHeight)));
        spatialFrame->vertices[3 * i + 0] = result.x;
        spatialFrame->vertices[3 * i + 1] = result.y;
        spatialFrame->vertices[3 * i + 2] = result.z;
    }
}

static void drawSpatialFrame(const Mesh spatialFrame)
{
    static const int frontFace[4][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
    static const int backFace[4][2] = {{4, 5}, {5, 6}, {6, 7}, {7, 4}};
    static const int connectingFaces[4][2] = {{0, 4}, {1, 7}, {2, 6}, {3, 5}};
    static const int (*faces[3])[2] = {frontFace, backFace, connectingFaces};
    Color near = NDC_SPACE() ? ANAKIWA : NEON_CARROT;
    Color ribs = ASPECT_CORRECT() ? MARINER : PALE_CANARY;
    Color far = PERSPECTIVE_CORRECT() ? ANAKIWA : HOPBUSH; //TODO: silly
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 4; j++)
        {
            const int (*edges)[2] = faces[i];
            const int edgeIndexA = edges[j][0] * 3;
            const int edgeIndexB = edges[j][1] * 3;
            const Vector3 startPosition = {spatialFrame.vertices[edgeIndexA + 0], spatialFrame.vertices[edgeIndexA + 1], spatialFrame.vertices[edgeIndexA + 2]};
            const Vector3 endPosition = {spatialFrame.vertices[edgeIndexB + 0], spatialFrame.vertices[edgeIndexB + 1], spatialFrame.vertices[edgeIndexB + 2]};
            const Color edgeColor = (i == 0) ? NEON_CARROT : (i == 1) ? EGGPLANT : HOPBUSH;
            // Color edgeColor = (i == 0) ? near : (i == 1) ? far : ribs;
            DrawLine3D(startPosition, endPosition, edgeColor);
        }
}

static void drawNearPlane(
    const Camera3D *main,
    const float aspect,
    const float near,
    const Model spatialFrameModel,
    const Mesh *mesh,
    const Texture2D meshTexture,
    const Texture2D perspectiveCorrectTexture,
    const float rotation)
{
    if (PERSPECTIVE_CORRECT() && TEXTURE_MODE())
    {
        spatialFrameModel.materials[0].maps[MATERIAL_MAP_ALBEDO].texture = perspectiveCorrectTexture;
        DrawModel(spatialFrameModel, MODEL_POS, 1.0f, WHITE);
        return;
    }
    perspectiveIncorrectCapture(main, aspect, near, mesh, meshTexture, rotation);
}

static Vector3
    aspectCorrectAndReflectNearPlane(const Vector3 hit, const Vector3 center, const Vector3 right, const Vector3 up, const float xAspect, const float yReflect)
{
    const Vector3 centerDistance = Vector3Subtract(hit, center);
    const float x = Vector3DotProduct(centerDistance, right);
    const float y = Vector3DotProduct(centerDistance, up);
    return Vector3Add(center, Vector3Add(Vector3Scale(right, x * xAspect), Vector3Scale(up, y * yReflect)));
}

static void
    drawNearPlanePoints(const Camera3D *main, const float aspect, const float near, const Model *nearPlanePointsModel, const Mesh *mesh, const float rotation)
{
    Mesh *nearPlanePointsMesh = &nearPlanePointsModel->meshes[0];
    const int capacity = mesh->triangleCount * 3;
    int writtenVertices = 0;
    Vector3 depth, right, up;
    basisVector(main, &depth, &right, &up);
    const Vector3 centerNearPlane = Vector3Add(main->position, Vector3Scale(depth, near));
    const float xAspect = Lerp(1.0f / aspect, 1.0f, aspectBlendFactor(0.0f));
    const float yReflect = Lerp(1.0f, -1.0f, reflectBlendFactor(0.0f));
    for (int i = 0; i < mesh->triangleCount; i++)
    {
        const int ia = mesh->indices ? mesh->indices[3 * i + 0] : 3 * i + 0;
        const int ib = mesh->indices ? mesh->indices[3 * i + 1] : 3 * i + 1;
        const int ic = mesh->indices ? mesh->indices[3 * i + 2] : 3 * i + 2;
        Vector3 a = (Vector3){mesh->vertices[3 * ia + 0], mesh->vertices[3 * ia + 1], mesh->vertices[3 * ia + 2]};
        Vector3 b = (Vector3){mesh->vertices[3 * ib + 0], mesh->vertices[3 * ib + 1], mesh->vertices[3 * ib + 2]};
        Vector3 c = (Vector3){mesh->vertices[3 * ic + 0], mesh->vertices[3 * ic + 1], mesh->vertices[3 * ic + 2]};
        a = translateRotateScale(0, a, MODEL_POS, MODEL_SCALE, rotation);
        b = translateRotateScale(0, b, MODEL_POS, MODEL_SCALE, rotation);
        c = translateRotateScale(0, c, MODEL_POS, MODEL_SCALE, rotation);
        const Vector3 *trianglePoints = (Vector3[]){a, b, c};
        //test if front facing or not
        if (Vector3DotProduct(Vector3Normalize(Vector3CrossProduct(Vector3Subtract(b, a), Vector3Subtract(c, a))), depth) > 0.0f) continue;

        const Vector3 hits[3] = {
            nearPlaneIntersection(main, near, a),
            nearPlaneIntersection(main, near, b),
            nearPlaneIntersection(main, near, c),
        };

        for (int j = 0; j < 3 && writtenVertices < capacity; ++j)
        {
            const Vector3 corrected = aspectCorrectAndReflectNearPlane(hits[j], centerNearPlane, right, up, xAspect, yReflect);
            DrawLine3D(trianglePoints[j], corrected, (Color){RED_DAMASK.r, RED_DAMASK.g, RED_DAMASK.b, 20});
            nearPlanePointsMesh->vertices[3 * writtenVertices + 0] = corrected.x;
            nearPlanePointsMesh->vertices[3 * writtenVertices + 1] = corrected.y;
            nearPlanePointsMesh->vertices[3 * writtenVertices + 2] = corrected.z;
            writtenVertices++;
        }
    }

    nearPlanePointsMesh->vertexCount = writtenVertices;
    rlSetPointSize(3.0f);
    DrawModelPoints(*nearPlanePointsModel, MODEL_POS, 1.0f, LILAC);
}

static void
    perspectiveIncorrectCapture(const Camera3D *main, const float aspect, const float near, const Mesh *mesh, const Texture2D meshTexture, const float rotation)
{
    Vector3 depth, right, up;
    basisVector(main, &depth, &right, &up);
    const Vector3 centerNearPlane = Vector3Add(main->position, Vector3Scale(depth, near));
    const float xAspect = Lerp(1.0f / aspect, 1.0f, aspectBlendFactor(0.0f));
    const float yReflect = Lerp(1.0f, -1.0f, reflectBlendFactor(0.0f));
    rlColor4ub(WHITE.r, WHITE.g, WHITE.b, WHITE.a);
    if (TEXTURE_MODE())
        rlEnableTexture(meshTexture.id);
    else
        rlDisableTexture();

    if (!TEXTURE_MODE() && !COLOR_MODE())
    {
        rlEnableWireMode();
        rlColor4ub(MARINER.r, MARINER.g, MARINER.b, MARINER.a);
    }
    rlBegin(RL_TRIANGLES);

    for (int i = 0; i < mesh->triangleCount; i++)
    {
        const int ia = mesh->indices ? mesh->indices[3 * i + 0] : 3 * i + 0;
        const int ib = mesh->indices ? mesh->indices[3 * i + 1] : 3 * i + 1;
        const int ic = mesh->indices ? mesh->indices[3 * i + 2] : 3 * i + 2;

        Vector3 vertexA = (Vector3){mesh->vertices[3 * ia + 0], mesh->vertices[3 * ia + 1], mesh->vertices[3 * ia + 2]};
        Vector3 vertexB = (Vector3){mesh->vertices[3 * ib + 0], mesh->vertices[3 * ib + 1], mesh->vertices[3 * ib + 2]};
        Vector3 vertexC = (Vector3){mesh->vertices[3 * ic + 0], mesh->vertices[3 * ic + 1], mesh->vertices[3 * ic + 2]};

        vertexA = translateRotateScale(0, vertexA, MODEL_POS, MODEL_SCALE, rotation);
        vertexB = translateRotateScale(0, vertexB, MODEL_POS, MODEL_SCALE, rotation);
        vertexC = translateRotateScale(0, vertexC, MODEL_POS, MODEL_SCALE, rotation);

        vertexA = aspectCorrectAndReflectNearPlane(nearPlaneIntersection(main, near, vertexA), centerNearPlane, right, up, xAspect, yReflect);
        vertexB = aspectCorrectAndReflectNearPlane(nearPlaneIntersection(main, near, vertexB), centerNearPlane, right, up, xAspect, yReflect);
        vertexC = aspectCorrectAndReflectNearPlane(nearPlaneIntersection(main, near, vertexC), centerNearPlane, right, up, xAspect, yReflect);

        if (COLOR_MODE()) rlColor4ub(mesh->colors[4 * ia + 0], mesh->colors[4 * ia + 1], mesh->colors[4 * ia + 2], mesh->colors[4 * ia + 3]);
        if (TEXTURE_MODE()) rlTexCoord2f(mesh->texcoords[2 * ia + 0], mesh->texcoords[2 * ia + 1]);
        rlVertex3f(vertexA.x, vertexA.y, vertexA.z);

        const int i2 = NDC_SPACE() && REFLECT_Y() ? ic : ib;
        const Vector3 vertex2 = NDC_SPACE() && REFLECT_Y() ? vertexC : vertexB;

        if (COLOR_MODE()) rlColor4ub(mesh->colors[4 * i2 + 0], mesh->colors[4 * i2 + 1], mesh->colors[4 * i2 + 2], mesh->colors[4 * i2 + 3]);
        if (TEXTURE_MODE()) rlTexCoord2f(mesh->texcoords[2 * i2 + 0], mesh->texcoords[2 * i2 + 1]);
        rlVertex3f(vertex2.x, vertex2.y, vertex2.z);

        const int i3 = NDC_SPACE() && REFLECT_Y() ? ib : ic;
        const Vector3 vertex3 = NDC_SPACE() && REFLECT_Y() ? vertexB : vertexC;

        if (COLOR_MODE()) rlColor4ub(mesh->colors[4 * i3 + 0], mesh->colors[4 * i3 + 1], mesh->colors[4 * i3 + 2], mesh->colors[4 * i3 + 3]);
        if (TEXTURE_MODE()) rlTexCoord2f(mesh->texcoords[2 * i3 + 0], mesh->texcoords[2 * i3 + 1]);
        rlVertex3f(vertex3.x, vertex3.y, vertex3.z);
    }

    rlEnd();
    rlDisableTexture();
    rlDisableWireMode();
}

static void
    perspectiveCorrectCapture(const Camera3D *main, const Model *model, const Texture2D meshTexture, Texture2D *perspectiveCorrectTexture, const float rotation)
{
    unsigned char *cacheVertexColors = model->meshes[0].colors;
    if (TEXTURE_MODE() && !COLOR_MODE()) model->meshes[0].colors = NULL;
    ClearBackground(BLACK);
    BeginMode3D(*main);
    //TODO: this is still needded to avoid the model in space and model proejected to near plane texture divergence... ugly though
    const Texture2D previousTexture = model->materials[0].maps[MATERIAL_MAP_ALBEDO].texture;
    model->materials[0].maps[MATERIAL_MAP_ALBEDO].texture = meshTexture;
    DrawModelEx(*model, MODEL_POS, (Vector3){0, 1, 0}, RAD2DEG * rotation, MODEL_SCALE, WHITE);
    model->materials[0].maps[MATERIAL_MAP_ALBEDO].texture = previousTexture;
    EndMode3D();
    Image rgba = LoadImageFromScreen();
    ImageFormat(&rgba, PIXELFORMAT_UNCOMPRESSED_R8G8B8A8);
    if (TEXTURE_MODE() && !COLOR_MODE()) model->meshes[0].colors = cacheVertexColors;
    ClearBackground(BLACK);
    BeginMode3D(*main);
    const Texture2D cacheTexture = model->materials[0].maps[MATERIAL_MAP_ALBEDO].texture;
    const Color cacheMaterialColor = model->materials[0].maps[MATERIAL_MAP_ALBEDO].color;
    model->materials[0].maps[MATERIAL_MAP_ALBEDO].texture = (Texture2D){0};
    model->materials[0].maps[MATERIAL_MAP_ALBEDO].color = WHITE;
    DrawModelEx(*model, MODEL_POS, (Vector3){0, 1, 0}, RAD2DEG * rotation, MODEL_SCALE, WHITE);
    model->materials[0].maps[MATERIAL_MAP_ALBEDO].texture = cacheTexture;
    model->materials[0].maps[MATERIAL_MAP_ALBEDO].color = cacheMaterialColor;
    EndMode3D();
    const Image mask = LoadImageFromScreen();
    alphaMaskPunchOut(&rgba, &mask, 1);
    ImageFlipVertical(&rgba);                                 //TODO: still a lot of pedagogy needs to be reorganized.... fix emphasis
    if (NDC_SPACE() && REFLECT_Y()) ImageFlipVertical(&rgba); //TODO: FLIP AGAIN... it works visually i think, but not clear enough, and also hacked/ugly
    if (perspectiveCorrectTexture->id != 0)
        UpdateTexture(*perspectiveCorrectTexture, rgba.data);
    else
        *perspectiveCorrectTexture = LoadTextureFromImage(rgba);

    UnloadImage(mask);
    UnloadImage(rgba);
}

static void alphaMaskPunchOut(Image *rgba, const Image *mask, const unsigned char threshold)
{
    Image maskGray = ImageCopy(*mask);
    ImageFormat(&maskGray, PIXELFORMAT_UNCOMPRESSED_GRAYSCALE);
    ImageFormat(rgba, PIXELFORMAT_UNCOMPRESSED_R8G8B8A8);
    const unsigned char *maskGrayScale = maskGray.data;
    Color *colors = rgba->data;
    const int pixelCount = rgba->width * rgba->height;
    for (size_t i = 0; i < pixelCount; ++i) colors[i].a = (maskGrayScale[i] > threshold) ? 255 : 0;
    UnloadImage(maskGray);
}

static float spaceBlendFactor(const float deltaTime)
{
    static float blend = 0.0f;
    if (deltaTime > 0.0f) blend = Clamp(blend + (NDC_SPACE() ? 1.0f : -1.0f) * BLEND_SCALAR * deltaTime, 0.0f, 1.0f);
    return blend;
}

static float aspectBlendFactor(const float deltaTime)
{
    static float blend = 0.0f;
    if (deltaTime > 0.0f) blend = Clamp(blend + (ASPECT_CORRECT() ? 1.0f : -1.0f) * BLEND_SCALAR * deltaTime, 0.0f, 1.0f);
    return blend;
}

static float reflectBlendFactor(const float deltaTime)
{
    static float blend = 0.0f;
    if (deltaTime > 0.0f)
    {
        const float target = (NDC_SPACE() && REFLECT_Y()) ? 1.0f : 0.0f;
        const float dir = (blend < target) ? 1.0f : (blend > target) ? -1.0f : 0.0f;
        blend = Clamp(blend + dir * BLEND_SCALAR * deltaTime, 0.0f, 1.0f);
    }
    return blend;
}

static void fillVertexColors(Mesh *mesh)
{
    if (!mesh->colors) mesh->colors = (unsigned char *)RL_CALLOC(mesh->vertexCount, sizeof(Color));
    Color *colors = (Color *)mesh->colors;
    const Vector3 *vertices = (Vector3 *)mesh->vertices;
    const BoundingBox bounds = GetMeshBoundingBox(*mesh);
    for (int i = 0; i < mesh->vertexCount; ++i)
    {
        const Vector3 vertex = vertices[i];
        float nx = (vertex.x - 0.5f * (bounds.min.x + bounds.max.x)) / (0.5f * (bounds.max.x - bounds.min.x));
        float ny = (vertex.y - 0.5f * (bounds.min.y + bounds.max.y)) / (0.5f * (bounds.max.y - bounds.min.y));
        float nz = (vertex.z - 0.5f * (bounds.min.z + bounds.max.z)) / (0.5f * (bounds.max.z - bounds.min.z));
        const float len = sqrtf(nx * nx + ny * ny + nz * nz);
        nx /= len;
        ny /= len;
        nz /= len;
        colors[i] =
            (Color){(unsigned char)lrintf(127.5f * (nx + 1.f)), (unsigned char)lrintf(127.5f * (ny + 1.f)), (unsigned char)lrintf(127.5f * (nz + 1.f)), 255};
    }
}
