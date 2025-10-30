#include "raylib.h"
#include "raymath.h"
#include "rlgl.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

//TODO:
// 1. add proper clipping to the target meshes to show how it works with the spaces
//   - move mesh out of clip planes or allow moving the main camera's target away from the meshes
// 2. improve didactic annotations (ideally with spatial labeling rather than simple flat screen overlay)
// 3. improve code didactic, code should read in order of fixed function staging
// 4. add scripted toggling/navigation of ordered fixed function staging visualization
// 5. add some sort of ghosting effect between fixed function stages, to emphasize previous stage perhaps)
// 6. OPTIONALLY improve toggling and space navigation

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

static const int FONT_SIZE = 20;
static const int WIDTH = 800;
static const int HEIGHT = 450;
static const float ANGULAR_VELOCITY = 1.25f;
static const float FOVY_PERSPECTIVE = 60.0f;
static const float BLEND_SCALAR = 5.0f;
static const Vector3 Y = {0.0f, 1.0f, 0.0f};
static const Vector3 MODEL_POS = {0.0f, 0.0f, 0.0f};
static const Vector3 MODEL_SCALE = {1.0f, 1.0f, 1.0f};
static const Vector3 MAIN_POS = {0.0f, 0.0f, 2.0f};
static const Vector3 JUGEMU_POS_ISO = {3.0f, 1.0f, 3.0f};

typedef unsigned short Triangle[3];
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

static void orbitSpace(Camera3D *jugemu, float dt);
static void basisVector(const Camera3D *main, Vector3 *depthOut, Vector3 *rightOut, Vector3 *upOut);
static void worldToNDCSpace(const Camera3D *main, float aspect, float near, float far, const Model *world, const Model *ndc, float rotation);
static void drawModelFilled(const Model *model, Texture2D texture, float rotation);
static void drawModelWiresAndPoints(const Model *model, float rotation);
static void drawNearPlanePoints(const Camera3D *main, float aspect, float near, const Model *nearPlanePointsModel, const Mesh *mesh, float rotation);
static void updateSpatialFrame(const Camera3D *main, float aspect, float near, float far, const Mesh *spatialFrame);
static void drawSpatialFrame(Mesh spatialFrame);
static void perspectiveIncorrectCapture(const Camera3D *main, float aspect, float near, const Mesh *mesh, Texture2D meshTexture, float rotation);
static void perspectiveCorrectCapture(const Camera3D *main, const Model *model, Texture2D meshTexture, Texture2D *perspectiveCorrectTexture, float rotation);
static void alphaMaskPunchOut(Image *rgba, const Image *mask, unsigned char threshold);
static void fillVertexColors(Mesh *mesh);
static float spaceBlendFactor(float dt);
static float aspectBlendFactor(float dt);
static float reflectBlendFactor(float dt);

int main(void)
{
    InitWindow(WIDTH, HEIGHT, "fixed function didactic");
    SetTargetFPS(60);
    Camera3D main = {0};
    main.position = MAIN_POS;
    main.target = MODEL_POS;
    main.up = Y;
    main.fovy = FOVY_PERSPECTIVE;
    main.projection = CAMERA_PERSPECTIVE;
    float aspect = (float)GetScreenWidth() / (float)GetScreenHeight();
    Camera3D jugemu = (Camera3D){0};
    jugemu.position = JUGEMU_POS_ISO;
    jugemu.target = MODEL_POS;
    jugemu.up = Y;
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

    if (!worldModel.meshes[0].indices)
    {
        worldModel.meshes[0].indices = RL_CALLOC(worldModel.meshes[0].vertexCount, sizeof(unsigned short));
        for (int i = 0; i < worldModel.meshes[0].vertexCount; i++) worldModel.meshes[0].indices[i] = (unsigned short)i;
        worldModel.meshes[0].triangleCount = worldModel.meshes[0].vertexCount / 3;
    }
    fillVertexColors(&worldModel.meshes[0]);

    Mesh ndcMesh = (Mesh){0};
    ndcMesh.vertexCount = worldModel.meshes[0].vertexCount;
    ndcMesh.triangleCount = worldModel.meshes[0].triangleCount;
    ndcMesh.vertices = RL_CALLOC(ndcMesh.vertexCount, sizeof(Vector3));
    ndcMesh.texcoords = RL_CALLOC(ndcMesh.vertexCount, sizeof(Vector2));
    ndcMesh.indices = RL_CALLOC(ndcMesh.triangleCount, sizeof(unsigned short[3]));
    ndcMesh.colors = RL_CALLOC(ndcMesh.vertexCount, sizeof(Color));
    memcpy(ndcMesh.colors, worldModel.meshes[0].colors, ndcMesh.vertexCount * sizeof(Color));
    memcpy(ndcMesh.texcoords, worldModel.meshes[0].texcoords, ndcMesh.vertexCount * sizeof(Vector2));
    memcpy(ndcMesh.indices, worldModel.meshes[0].indices, ndcMesh.triangleCount * sizeof(unsigned short[3]));
    Model ndcModel = LoadModelFromMesh(ndcMesh);

    Mesh nearPlanePoints = (Mesh){0};
    nearPlanePoints.vertexCount = worldModel.meshes[0].triangleCount * 3;
    nearPlanePoints.vertices = RL_CALLOC(nearPlanePoints.vertexCount, sizeof(Vector3));
    Model nearPlanePointsModel = LoadModelFromMesh(nearPlanePoints);

    worldModel.materials[0].maps[MATERIAL_MAP_ALBEDO].texture = meshTexture;
    ndcModel.materials[0].maps[MATERIAL_MAP_ALBEDO].texture = meshTexture;

    Texture2D perspectiveCorrectTexture = {0};
    Mesh spatialFrame = GenMeshCube(1.0f, 1.0f, 1.0f);
    spatialFrame.colors = RL_CALLOC(spatialFrame.vertexCount, sizeof(Color));
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

        const float sBlend = spaceBlendFactor(GetFrameTime());
        aspectBlendFactor(GetFrameTime());
        reflectBlendFactor(GetFrameTime());

        if (!PAUSED()) meshRotation -= ANGULAR_VELOCITY * GetFrameTime();

        orbitSpace(&jugemu, GetFrameTime());

        worldToNDCSpace(&main, aspect, near, far, &worldModel, &ndcModel, meshRotation);

        for (int i = 0; i < ndcModel.meshes[0].vertexCount; ++i)
        {
            const Vector3 *worldVertices = (const Vector3 *)worldModel.meshes[0].vertices;
            Vector3 *ndcVertices = (Vector3 *)ndcModel.meshes[0].vertices;
            ndcVertices[i].x = Lerp(worldVertices[i].x, ndcVertices[i].x, sBlend);
            ndcVertices[i].y = Lerp(worldVertices[i].y, ndcVertices[i].y, sBlend);
            ndcVertices[i].z = Lerp(worldVertices[i].z, ndcVertices[i].z, sBlend);
        }
        Model *displayModel = &ndcModel;
        Mesh *displayMesh = &ndcModel.meshes[0];
        if (PERSPECTIVE_CORRECT() && TEXTURE_MODE()) perspectiveCorrectCapture(&main, displayModel, meshTexture, &perspectiveCorrectTexture, meshRotation);

        BeginDrawing();

        ClearBackground(BLACK);
        BeginMode3D(jugemu);
        Vector3 depth, right, up;
        basisVector(&main, &depth, &right, &up);
        DrawLine3D(main.position, Vector3Add(main.position, right), NEON_CARROT);
        DrawLine3D(main.position, Vector3Add(main.position, up), LILAC);
        DrawLine3D(main.position, Vector3Add(main.position, depth), MARINER);

        updateSpatialFrame(&main, aspect, near, far, &spatialFrame);
        drawSpatialFrame(spatialFrame);

        drawModelFilled(displayModel, meshTexture, meshRotation);
        drawModelWiresAndPoints(displayModel, meshRotation);

        drawNearPlanePoints(&main, aspect, near, &nearPlanePointsModel, displayMesh, meshRotation);
        if (PERSPECTIVE_CORRECT() && TEXTURE_MODE())
        {
            spatialFrameModel.materials[0].maps[MATERIAL_MAP_ALBEDO].texture = perspectiveCorrectTexture;
            DrawModel(spatialFrameModel, MODEL_POS, 1.0f, WHITE);
        }
        else
        {
            perspectiveIncorrectCapture(&main, aspect, near, displayMesh, meshTexture, meshRotation);
        }
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

static void orbitSpace(Camera3D *jugemu, const float dt)
{
    float radius = Vector3Length(jugemu->position);
    float azimuth = atan2f(jugemu->position.z, jugemu->position.x);
    const float horizontalRadius = sqrtf(jugemu->position.x * jugemu->position.x + jugemu->position.z * jugemu->position.z);
    float elevation = atan2f(jugemu->position.y, horizontalRadius);
    if (IsKeyDown(KEY_LEFT)) azimuth += 1.5f * dt;
    if (IsKeyDown(KEY_RIGHT)) azimuth -= 1.5f * dt;
    if (IsKeyDown(KEY_UP)) elevation += 1.0f * dt;
    if (IsKeyDown(KEY_DOWN)) elevation -= 1.0f * dt;
    if (IsKeyDown(KEY_W)) radius -= 2.0f * dt;
    if (IsKeyDown(KEY_S)) radius += 2.0f * dt;
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

static Vector3 translateRotateScale(const int inverse, const Vector3 coordinate, const Vector3 pos, const Vector3 scale, const float rotation)
{
    const Matrix matrix = MatrixMultiply(MatrixMultiply(MatrixScale(scale.x, scale.y, scale.z), MatrixRotateY(rotation)), MatrixTranslate(pos.x, pos.y, pos.z));
    const Matrix result = inverse ? MatrixInvert(matrix) : matrix;
    return Vector3Transform(coordinate, result);
}

static Vector3 intersect(const Camera3D *main, const float near, const Vector3 worldCoord)
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
    for (int i = 0; i < world->meshes[0].vertexCount; i++)
    {
        const Vector3 worldVertex = translateRotateScale(0, ((const Vector3 *)world->meshes[0].vertices)[i], MODEL_POS, MODEL_SCALE, rotation);
        const float signedDepth = Vector3DotProduct(Vector3Subtract(worldVertex, main->position), depth);
        const Vector3 intersectionCoord = intersect(main, near, worldVertex);
        const Vector3 clipPlaneVector = Vector3Subtract(intersectionCoord, centerNearPlane);
        const float xNDC = Vector3DotProduct(clipPlaneVector, right) / halfWNear;
        const float yNDC = Vector3DotProduct(clipPlaneVector, up) / halfHNear;
        const float zNDC = (far + near - 2.0f * far * near / signedDepth) / (far - near);
        const Vector3 scaledRight = Vector3Scale(right, xNDC * halfWNear);
        const Vector3 scaledUp = Vector3Scale(up, yNDC * halfHNear);
        const Vector3 scaledDepth = Vector3Scale(depth, zNDC * halfDepthNDC);
        const Vector3 offset = Vector3Add(Vector3Add(scaledRight, scaledUp), scaledDepth);
        const Vector3 scaledNDCCoord = Vector3Add(centerNDCCube, offset);
        ((Vector3 *)ndc->meshes[0].vertices)[i] = translateRotateScale(1, scaledNDCCoord, MODEL_POS, MODEL_SCALE, rotation);
    }
}

static void drawModelFilled(const Model *model, const Texture2D texture, const float rotation)
{
    if (!(COLOR_MODE() || TEXTURE_MODE())) return;
    Color *cacheColors = (Color *)model->meshes[0].colors;
    if (TEXTURE_MODE() && !COLOR_MODE()) model->meshes[0].colors = NULL;
    model->materials[0].maps[MATERIAL_MAP_ALBEDO].texture.id = TEXTURE_MODE() ? texture.id : 0;
    DrawModelEx(*model, MODEL_POS, Y, RAD2DEG * rotation, MODEL_SCALE, WHITE);
    model->materials[0].maps[MATERIAL_MAP_ALBEDO].texture.id = 0;
    model->meshes[0].colors = (unsigned char *)cacheColors;
}

static void drawModelWiresAndPoints(const Model *model, const float rotation)
{
    Color *cacheColors = (Color *)model->meshes[0].colors;
    const unsigned int cacheID = model->materials[0].maps[MATERIAL_MAP_ALBEDO].texture.id;
    model->meshes[0].colors = NULL;
    model->materials[0].maps[MATERIAL_MAP_ALBEDO].texture.id = 0;
    DrawModelWiresEx(*model, MODEL_POS, Y, RAD2DEG * rotation, MODEL_SCALE, MARINER);
    rlSetPointSize(4.0f);
    DrawModelPointsEx(*model, MODEL_POS, Y, RAD2DEG * rotation, MODEL_SCALE, LILAC);
    model->materials[0].maps[MATERIAL_MAP_ALBEDO].texture.id = cacheID;
    model->meshes[0].colors = (unsigned char *)cacheColors;
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
    const float halfDepthNdc = Lerp(halfHNear, 0.5f * (far - near), aspectBlendFactor(0.0f));
    const float halfDepth = Lerp(0.5f * (far - near), halfDepthNdc, spaceBlendFactor(0.0f));
    const float farHalfW = Lerp(halfWFar, halfWNear, spaceBlendFactor(0.0f));
    const float farHalfH = Lerp(halfHFar, halfHNear, spaceBlendFactor(0.0f));

    for (int i = 0; i < spatialFrame->vertexCount; ++i)
    {
        const Vector3 offset = Vector3Subtract(((Vector3 *)spatialFrame->vertices)[i], centerNear);
        const float xSign = Vector3DotProduct(offset, right) >= 0.0f ? 1.0f : -1.0f;
        const float ySign = Vector3DotProduct(offset, up) >= 0.0f ? 1.0f : -1.0f;
        const float farMask = Vector3DotProduct(offset, depth) > halfDepth ? 1.0f : 0.0f;
        const float finalHalfW = halfWNear + farMask * (farHalfW - halfWNear);
        const float finalHalfH = halfHNear + farMask * (farHalfH - halfHNear);
        const Vector3 center = Vector3Add(centerNear, Vector3Scale(depth, farMask * 2.0f * halfDepth));
        ((Vector3 *)spatialFrame->vertices)[i] = Vector3Add(center, Vector3Add(Vector3Scale(right, xSign * finalHalfW), Vector3Scale(up, ySign * finalHalfH)));
    }
}

static void drawSpatialFrame(const Mesh spatialFrame)
{
    static const int frontFace[4][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
    static const int backFace[4][2] = {{4, 5}, {5, 6}, {6, 7}, {7, 4}};
    static const int connectingFaces[4][2] = {{0, 4}, {1, 7}, {2, 6}, {3, 5}};
    static const int (*faces[3])[2] = {frontFace, backFace, connectingFaces};
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 4; j++)
        {
            const Vector3 startPosition = ((const Vector3 *)spatialFrame.vertices)[faces[i][j][0]];
            const Vector3 endPosition = ((const Vector3 *)spatialFrame.vertices)[faces[i][j][1]];
            DrawLine3D(startPosition, endPosition, i == 0 ? NEON_CARROT : i == 1 ? EGGPLANT : HOPBUSH);
        }
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
    Vector3 depth, right, up;
    basisVector(main, &depth, &right, &up);
    int nearPlaneVertexCount = 0;
    const int capacity = mesh->triangleCount * 3;
    Mesh *nearPlanePointsMesh = &nearPlanePointsModel->meshes[0];
    const Vector3 centerNearPlane = Vector3Add(main->position, Vector3Scale(depth, near));
    const float xAspect = Lerp(1.0f / aspect, 1.0f, aspectBlendFactor(0.0f));
    const float yReflect = Lerp(1.0f, -1.0f, reflectBlendFactor(0.0f));
    for (int i = 0; i < mesh->triangleCount; i++)
    {
        const Vector3 *vertices = (Vector3 *)mesh->vertices;
        const Triangle *triangles = (const Triangle *)mesh->indices;
        const Vector3 a = translateRotateScale(0, vertices[triangles[i][0]], MODEL_POS, MODEL_SCALE, rotation);
        const Vector3 b = translateRotateScale(0, vertices[triangles[i][1]], MODEL_POS, MODEL_SCALE, rotation);
        const Vector3 c = translateRotateScale(0, vertices[triangles[i][2]], MODEL_POS, MODEL_SCALE, rotation);
        // test if front facing or not (ugly one-liner -- comment out will ~double the rays)
        if (Vector3DotProduct(Vector3Normalize(Vector3CrossProduct(Vector3Subtract(b, a), Vector3Subtract(c, a))), depth) > 0.0f) continue;

        const Vector3 intersectionPoints[3] = {intersect(main, near, a), intersect(main, near, b), intersect(main, near, c)};
        for (int j = 0; j < 3 && nearPlaneVertexCount < capacity; ++j)
        {
            const Vector3 corrected = aspectCorrectAndReflectNearPlane(intersectionPoints[j], centerNearPlane, right, up, xAspect, yReflect);
            DrawLine3D((Vector3[]){a, b, c}[j], corrected, (Color){RED_DAMASK.r, RED_DAMASK.g, RED_DAMASK.b, 20});
            ((Vector3 *)nearPlanePointsMesh->vertices)[nearPlaneVertexCount] = corrected;
            nearPlaneVertexCount++;
        }
    }

    nearPlanePointsMesh->vertexCount = nearPlaneVertexCount;
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
        const Triangle *triangles = (const Triangle *)mesh->indices;
        const Vector3 *vertices = (const Vector3 *)mesh->vertices;
        const Color *colors = (const Color *)mesh->colors;
        const Vector2 *texcoords = (const Vector2 *)mesh->texcoords;

        Vector3 a = translateRotateScale(0, vertices[triangles[i][0]], MODEL_POS, MODEL_SCALE, rotation);
        Vector3 b = translateRotateScale(0, vertices[triangles[i][1]], MODEL_POS, MODEL_SCALE, rotation);
        Vector3 c = translateRotateScale(0, vertices[triangles[i][2]], MODEL_POS, MODEL_SCALE, rotation);

        a = aspectCorrectAndReflectNearPlane(intersect(main, near, a), centerNearPlane, right, up, xAspect, yReflect);
        b = aspectCorrectAndReflectNearPlane(intersect(main, near, b), centerNearPlane, right, up, xAspect, yReflect);
        c = aspectCorrectAndReflectNearPlane(intersect(main, near, c), centerNearPlane, right, up, xAspect, yReflect);

        if (COLOR_MODE()) rlColor4ub(colors[triangles[i][0]].r, colors[triangles[i][0]].g, colors[triangles[i][0]].b, colors[triangles[i][0]].a);
        if (TEXTURE_MODE()) rlTexCoord2f(texcoords[triangles[i][0]].x, texcoords[triangles[i][0]].y);
        rlVertex3f(a.x, a.y, a.z);
        // index winding to account for reflection toggle (will draw the inside of the geometry otherwise)
        const int secondIndex = NDC_SPACE() && REFLECT_Y() ? triangles[i][2] : triangles[i][1];
        const Vector3 secondVertex = NDC_SPACE() && REFLECT_Y() ? c : b;
        if (COLOR_MODE()) rlColor4ub(colors[secondIndex].r, colors[secondIndex].g, colors[secondIndex].b, colors[secondIndex].a);
        if (TEXTURE_MODE()) rlTexCoord2f(texcoords[secondIndex].x, texcoords[secondIndex].y);
        rlVertex3f(secondVertex.x, secondVertex.y, secondVertex.z);

        const int thirdIndex = NDC_SPACE() && REFLECT_Y() ? triangles[i][1] : triangles[i][2];
        const Vector3 thirdVertex = NDC_SPACE() && REFLECT_Y() ? b : c;
        if (COLOR_MODE()) rlColor4ub(colors[thirdIndex].r, colors[thirdIndex].g, colors[thirdIndex].b, colors[thirdIndex].a);
        if (TEXTURE_MODE()) rlTexCoord2f(texcoords[thirdIndex].x, texcoords[thirdIndex].y);
        rlVertex3f(thirdVertex.x, thirdVertex.y, thirdVertex.z);
    }

    rlEnd();
    rlDisableTexture();
    rlDisableWireMode();
}

static void
    perspectiveCorrectCapture(const Camera3D *main, const Model *model, const Texture2D meshTexture, Texture2D *perspectiveCorrectTexture, const float rotation)
{
    unsigned char *cacheColors = model->meshes[0].colors;
    if (TEXTURE_MODE() && !COLOR_MODE()) model->meshes[0].colors = NULL;

    ClearBackground(BLACK);
    BeginMode3D(*main);
    const Texture2D previousTexture = model->materials[0].maps[MATERIAL_MAP_ALBEDO].texture;
    model->materials[0].maps[MATERIAL_MAP_ALBEDO].texture = meshTexture;
    DrawModelEx(*model, MODEL_POS, Y, RAD2DEG * rotation, MODEL_SCALE, WHITE);
    model->materials[0].maps[MATERIAL_MAP_ALBEDO].texture = previousTexture;
    EndMode3D();

    Image rgba = LoadImageFromScreen();
    ImageFormat(&rgba, PIXELFORMAT_UNCOMPRESSED_R8G8B8A8);
    model->meshes[0].colors = cacheColors;

    ClearBackground(BLACK);
    BeginMode3D(*main);
    const Texture2D cacheTexture = model->materials[0].maps[MATERIAL_MAP_ALBEDO].texture;
    const Color cacheMaterialColor = model->materials[0].maps[MATERIAL_MAP_ALBEDO].color;
    model->materials[0].maps[MATERIAL_MAP_ALBEDO].texture = (Texture2D){0};
    model->materials[0].maps[MATERIAL_MAP_ALBEDO].color = WHITE;
    DrawModelEx(*model, MODEL_POS, Y, RAD2DEG * rotation, MODEL_SCALE, WHITE);
    model->materials[0].maps[MATERIAL_MAP_ALBEDO].texture = cacheTexture;
    model->materials[0].maps[MATERIAL_MAP_ALBEDO].color = cacheMaterialColor;
    EndMode3D();

    const Image mask = LoadImageFromScreen();
    alphaMaskPunchOut(&rgba, &mask, 1);
    ImageFlipVertical(&rgba);
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
    Image maskCopy = ImageCopy(*mask);
    ImageFormat(&maskCopy, PIXELFORMAT_UNCOMPRESSED_GRAYSCALE);
    ImageFormat(rgba, PIXELFORMAT_UNCOMPRESSED_R8G8B8A8);
    const unsigned char *maskGrayScale = maskCopy.data;
    Color *colors = rgba->data;
    const int pixelCount = rgba->width * rgba->height;
    for (size_t i = 0; i < pixelCount; ++i) colors[i].a = maskGrayScale[i] > threshold ? 255 : 0;
    UnloadImage(maskCopy);
}

static float spaceBlendFactor(const float dt)
{
    static float blend = 0.0f;
    if (dt > 0.0f) blend = Clamp(blend + (NDC_SPACE() ? 1.0f : -1.0f) * BLEND_SCALAR * dt, 0.0f, 1.0f);
    return blend;
}

static float aspectBlendFactor(const float dt)
{
    static float blend = 0.0f;
    if (dt > 0.0f) blend = Clamp(blend + (ASPECT_CORRECT() ? 1.0f : -1.0f) * BLEND_SCALAR * dt, 0.0f, 1.0f);
    return blend;
}

static float reflectBlendFactor(const float dt)
{
    static float blend = 0.0f;
    if (dt > 0.0f)
    {
        const float target = NDC_SPACE() && REFLECT_Y() ? 1.0f : 0.0f;
        const float direction = blend < target ? 1.0f : blend > target ? -1.0f : 0.0f;
        blend = Clamp(blend + direction * BLEND_SCALAR * dt, 0.0f, 1.0f);
    }
    return blend;
}

static void fillVertexColors(Mesh *mesh)
{
    if (!mesh->colors) mesh->colors = RL_CALLOC(mesh->vertexCount, sizeof(Color));
    Color *colors = (Color *)mesh->colors;
    const Vector3 *vertices = (Vector3 *)mesh->vertices;
    const BoundingBox bounds = GetMeshBoundingBox(*mesh);

    for (int i = 0; i < mesh->vertexCount; ++i)
    {
        const Vector3 vertex = vertices[i];
        const float nx = (vertex.x - 0.5f * (bounds.min.x + bounds.max.x)) / (0.5f * (bounds.max.x - bounds.min.x));
        const float ny = (vertex.y - 0.5f * (bounds.min.y + bounds.max.y)) / (0.5f * (bounds.max.y - bounds.min.y));
        const float nz = (vertex.z - 0.5f * (bounds.min.z + bounds.max.z)) / (0.5f * (bounds.max.z - bounds.min.z));
        const float len = sqrtf(nx * nx + ny * ny + nz * nz);
        colors[i] = (Color){lrintf(127.5f * (nx / len + 1.f)), lrintf(127.5f * (ny / len + 1.f)), lrintf(127.5f * (nz / len + 1.f)), 255};
    }
}
