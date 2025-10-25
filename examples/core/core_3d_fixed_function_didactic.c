#include "raylib.h"
#include "raymath.h"
#include "rlgl.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

static const int WIDTH = 640;
static const int HEIGHT = 480;
static const float ANGULAR_VELOCITY = 1.25f;
static const float FOVY_PERSPECTIVE = 60.0f;
static const Vector3 MODEL_POS = {0.0f, 0.0f, 0.0f};
static const Vector3 MODEL_SCALE = {1.0f, 1.0f, 1.0f};
static const Vector3 MAIN_POS = {0.0f, 0.0f, 2.0f};
static const Vector3 JUGEMU_POS_ISO = {3.0f, 1.0f, 3.0f};

enum
{
    FLAG_NDC_OVERLAY = 1u << 0,         // N
    FLAG_ASPECT = 1u << 1,              // Q
    FLAG_PERSPECTIVE_CORRECT = 1u << 2, // P
    FLAG_PAUSE = 1u << 3,               // F
    FLAG_COLOR_MODE = 1u << 4,          // C
    FLAG_TEXTURE_MODE = 1u << 5         // T
};
static unsigned int gFlags = FLAG_ASPECT | FLAG_COLOR_MODE;

#define NDC_SPACE() ((gFlags & FLAG_NDC_OVERLAY) != 0)
#define ANISOTROPIC() ((gFlags & FLAG_ASPECT) != 0)
#define PERSPECTIVE_CORRECT() ((gFlags & FLAG_PERSPECTIVE_CORRECT) != 0)
#define PAUSED() ((gFlags & FLAG_PAUSE) != 0)
#define COLOR_MODE() ((gFlags & FLAG_COLOR_MODE) != 0)
#define TEXTURE_MODE() ((gFlags & FLAG_TEXTURE_MODE) != 0)

static void drawObservedAxes(Camera3D *main);
static void drawFrustum(Camera3D *main, float aspect, float near, float far);
static void drawNearPlane(Camera3D *main, float aspect, float near, Mesh *mesh, Texture2D texture, float rotation);

static void drawModelFilled(Model *model, Texture2D texture, float rotation);
static void drawWiresAndPoints(Model *model, float rotation);

static int nearPlanePointsCapacity = 0;
static void drawNearPlanePoints(Camera3D *main, float near, Model *nearPlanePointsModel, Mesh *mesh, float rotation, bool flipY);

static void worldToNDCSpace(Camera3D *main, float aspect, float near, float far, Model *world, Model *ndc, float rotation);
static void worldToNDCMesh(Camera3D *main, float aspect, float near, float far, Mesh *world, Mesh *ndc, float rotation);

static Mesh nearPlaneQuad = {0};
static Model nearPlaneQuadModel = {0};
static void buildNearPlaneQuad(void);
static void updateNearPlaneQuad(Camera3D *main, float aspect, float near);

static void perspectiveIncorrectTextureCapture(Camera3D *main, float aspect, float near, Mesh *mesh, Texture2D texture, float rotation);
static void flatColor(Camera3D *main, float aspect, float near, Mesh *mesh, float rotation);
static void applyBarycentricPalette(Mesh *mesh);

static Texture2D perspectiveCorrectTexture = {0};
static void perspectiveCorrectTextureCapture(Camera3D *main, Model *model, Texture2D texture, float rotation);
static void punchOutAlphaFromImage(Image *rgba, const Image *mask, unsigned char threshold);
static void moveJugemuOrbital(Camera3D *jugemu, float deltaTime);

int main(void)
{
    InitWindow(WIDTH, HEIGHT, "fixed function");
    SetTargetFPS(60);

    Camera3D main = {0};
    main.position = MAIN_POS;
    main.target = (Vector3){0, 0, 0};
    main.up = (Vector3){0, 1, 0};
    main.fovy = FOVY_PERSPECTIVE;
    main.projection = CAMERA_PERSPECTIVE;

    float near = 1.0f;
    float far = 3.0f;
    float aspect = (float)GetScreenWidth() / (float)GetScreenHeight();

    Camera3D jugemu = (Camera3D){0};
    jugemu.position = JUGEMU_POS_ISO;
    jugemu.target = (Vector3){0, 0, 0};
    jugemu.up = (Vector3){0, 1, 0};
    jugemu.fovy = FOVY_PERSPECTIVE;
    jugemu.projection = CAMERA_PERSPECTIVE;

    float meshRotation = 0.0f;
    Mesh cubeMesh = GenMeshCube(1.0f, 1.0f, 1.0f);
    Model worldModel = LoadModelFromMesh(cubeMesh);

    Mesh ndcMesh = (Mesh){0};
    ndcMesh.vertexCount = worldModel.meshes[0].vertexCount;
    ndcMesh.triangleCount = worldModel.meshes[0].triangleCount;
    ndcMesh.vertices = RL_CALLOC(ndcMesh.vertexCount * 3, sizeof(float));
    ndcMesh.texcoords = RL_CALLOC(ndcMesh.vertexCount * 2, sizeof(float));
    ndcMesh.indices = RL_CALLOC(ndcMesh.triangleCount * 3, sizeof(unsigned short));

    if (worldModel.meshes[0].texcoords)
        memcpy(ndcMesh.texcoords, worldModel.meshes[0].texcoords, sizeof(float) * 2 * ndcMesh.vertexCount);

    if (worldModel.meshes[0].indices)
        memcpy(ndcMesh.indices, worldModel.meshes[0].indices, sizeof(unsigned short) * 3 * ndcMesh.triangleCount);

    Model ndcModel = LoadModelFromMesh(ndcMesh);

    Mesh nearPlanePoints = (Mesh){0};
    nearPlanePoints.vertexCount = 3 * worldModel.meshes[0].triangleCount;
    nearPlanePointsCapacity = nearPlanePoints.vertexCount;
    nearPlanePoints.vertices = RL_CALLOC(nearPlanePoints.vertexCount * 3, sizeof(float));
    Model nearPlanePointsModel = LoadModelFromMesh(nearPlanePoints);

    Image checkeredImage = GenImageChecked(16, 16, 4, 4, BLACK, WHITE);
    Texture2D texture = LoadTextureFromImage(checkeredImage);
    UnloadImage(checkeredImage);

    worldModel.materials[0].maps[MATERIAL_MAP_ALBEDO].texture = texture;
    applyBarycentricPalette(&worldModel.meshes[0]);
    ndcModel.materials[0].maps[MATERIAL_MAP_ALBEDO].texture = texture;
    applyBarycentricPalette(&ndcModel.meshes[0]);

    while (!WindowShouldClose())
    {
        aspect = (float)GetScreenWidth() / (float)GetScreenHeight();
        if (IsKeyPressed(KEY_N))
            gFlags ^= FLAG_NDC_OVERLAY;
        if (IsKeyPressed(KEY_Q))
            gFlags ^= FLAG_ASPECT;
        if (IsKeyPressed(KEY_P))
            gFlags ^= FLAG_PERSPECTIVE_CORRECT;
        if (IsKeyPressed(KEY_F))
            gFlags ^= FLAG_PAUSE;
        if (IsKeyPressed(KEY_C))
            gFlags ^= FLAG_COLOR_MODE;
        if (IsKeyPressed(KEY_T))
            gFlags ^= FLAG_TEXTURE_MODE;

        moveJugemuOrbital(&jugemu, GetFrameTime());
        if (!PAUSED())
            meshRotation -= ANGULAR_VELOCITY * GetFrameTime();

        BeginDrawing();
        ClearBackground(BLACK);
        if (TEXTURE_MODE() && PERSPECTIVE_CORRECT())
        {
            if (NDC_SPACE())
            {
                worldToNDCMesh(&main, aspect, near, far, &worldModel.meshes[0], &ndcModel.meshes[0], meshRotation);
                perspectiveCorrectTextureCapture(&main, &ndcModel, texture, meshRotation);
            }
            else
            {
                perspectiveCorrectTextureCapture(&main, &worldModel, texture, meshRotation);
            }
            ClearBackground(BLACK);
        }

        BeginMode3D(jugemu);
        drawObservedAxes(&main);
        if (NDC_SPACE())
        {
            worldToNDCSpace(&main, aspect, near, far, &worldModel, &ndcModel, meshRotation);
            drawModelFilled(&ndcModel, texture, meshRotation);
            drawWiresAndPoints(&ndcModel, meshRotation);
            drawNearPlanePoints(&main, near, &nearPlanePointsModel, &ndcModel.meshes[0], meshRotation, true);
            drawNearPlane(&main, aspect, near, &ndcModel.meshes[0], texture, meshRotation);
        }
        else
        {
            drawModelFilled(&worldModel, texture, meshRotation);
            drawWiresAndPoints(&worldModel, meshRotation);
            drawFrustum(&main, ANISOTROPIC() ? aspect : 1.0f, near, far);
            drawNearPlanePoints(&main, near, &nearPlanePointsModel, &worldModel.meshes[0], meshRotation, false);
            drawNearPlane(&main, aspect, near, &worldModel.meshes[0], texture, meshRotation);
        }
        EndMode3D();
        DrawText(
            "E: NDC | Q: aspect | P: perspective divide | C: colors | T: texture | arrowkeys/WASD: move camera | Space: reset camera", 12, 14, 10, RAYWHITE);

        int lineY = 32;
        int lineStep = 12;
        DrawText(TextFormat("SPACE:   %s", NDC_SPACE() ? "NDC" : "WORLD"), 12, lineY, 10, RAYWHITE);
        lineY += lineStep;
        DrawText(TextFormat("ASPECT:  %s", ANISOTROPIC() ? "ANISOTROPIC" : "ISOTROPIC"), 12, lineY, 10, RAYWHITE);
        lineY += lineStep;
        DrawText(TextFormat("TEXTURE: %s", TEXTURE_MODE() ? "ON" : "OFF"), 12, lineY, 10, RAYWHITE);
        lineY += lineStep;
        DrawText(TextFormat("PERSPECTIVE:  %s", PERSPECTIVE_CORRECT() ? "CORRECT" : "INCORRECT"), 12, lineY, 10, RAYWHITE);
        lineY += lineStep;
        DrawText(TextFormat("COLORS:  %s", COLOR_MODE() ? "ON" : "OFF"), 12, lineY, 10, RAYWHITE);
        EndDrawing();
    }

    UnloadModel(worldModel);
    UnloadModel(ndcModel);
    UnloadModel(nearPlanePointsModel);
    UnloadTexture(texture);
    CloseWindow();
    return 0;
}

static void basisVector(Camera3D *main, Vector3 *losOut, Vector3 *rightOut, Vector3 *upOut)
{
    Vector3 los = Vector3Normalize(Vector3Subtract(main->target, main->position));
    Vector3 right = Vector3Normalize(Vector3CrossProduct(los, main->up));
    Vector3 up = Vector3Normalize(Vector3CrossProduct(right, los));
    *losOut = los;
    *rightOut = right;
    *upOut = up;
}

static void drawObservedAxes(Camera3D *main)
{
    Vector3 los, right, up;
    basisVector(main, &los, &right, &up);
    DrawLine3D(main->position, Vector3Add(main->position, right), PURPLE);
    DrawLine3D(main->position, Vector3Add(main->position, up), DARKGREEN);
    DrawLine3D(main->position, Vector3Add(main->position, los), SKYBLUE);
}

static void drawFrustum(Camera3D *main, float aspect, float near, float far)
{
    //TODO: there is likely an easier more elegant way to do this, like a mesh GEN or something and modify not just drawing the lines
    Vector3 los, right, up;
    basisVector(main, &los, &right, &up);
    Vector3 centerNear = Vector3Add(main->position, Vector3Scale(los, near));
    Vector3 centerFar = Vector3Add(main->position, Vector3Scale(los, far));
    float halfFovy = DEG2RAD * main->fovy * 0.5f;
    float halfHNear = near * tanf(halfFovy);
    float halfWNear = halfHNear * aspect;
    float halfHFar = far * tanf(halfFovy);
    float halfWFar = halfHFar * aspect;
    Vector3 nearTL = Vector3Add(Vector3Add(centerNear, Vector3Scale(up, halfHNear)), Vector3Scale(right, -halfWNear));
    Vector3 nearTR = Vector3Add(Vector3Add(centerNear, Vector3Scale(up, halfHNear)), Vector3Scale(right, halfWNear));
    Vector3 nearBR = Vector3Add(Vector3Add(centerNear, Vector3Scale(up, -halfHNear)), Vector3Scale(right, halfWNear));
    Vector3 nearBL = Vector3Add(Vector3Add(centerNear, Vector3Scale(up, -halfHNear)), Vector3Scale(right, -halfWNear));
    Vector3 farTL = Vector3Add(Vector3Add(centerFar, Vector3Scale(up, halfHFar)), Vector3Scale(right, -halfWFar));
    Vector3 farTR = Vector3Add(Vector3Add(centerFar, Vector3Scale(up, halfHFar)), Vector3Scale(right, halfWFar));
    Vector3 farBR = Vector3Add(Vector3Add(centerFar, Vector3Scale(up, -halfHFar)), Vector3Scale(right, halfWFar));
    Vector3 farBL = Vector3Add(Vector3Add(centerFar, Vector3Scale(up, -halfHFar)), Vector3Scale(right, -halfWFar));
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

static void drawNDCCube(Camera3D *main, Vector3 centerNear, float halfWNear, float halfHNear, float halfDepthNDCCube)
{
    Vector3 los, right, up;
    basisVector(main, &los, &right, &up);
    Vector3 corners[8];
    int writeIndex = 0;

    for (int zSign = -1; zSign <= 1; zSign += 2)
        for (int ySign = -1; ySign <= 1; ySign += 2)
            for (int xSign = -1; xSign <= 1; xSign += 2)
            {
                corners[writeIndex++] = Vector3Add(
                    centerNear,
                    Vector3Add(
                        Vector3Scale(right, xSign * halfWNear), Vector3Add(Vector3Scale(up, ySign * halfHNear), Vector3Scale(los, zSign * halfDepthNDCCube))));
            }

    Vector3 nearTL = corners[3], nearTR = corners[2], nearBR = corners[0], nearBL = corners[1];
    Vector3 farTL = corners[7], farTR = corners[6], farBR = corners[4], farBL = corners[5];

    DrawLine3D(nearTL, nearTR, SKYBLUE);
    DrawLine3D(nearTR, nearBR, SKYBLUE);
    DrawLine3D(nearBR, nearBL, SKYBLUE);
    DrawLine3D(nearBL, nearTL, SKYBLUE);
    DrawLine3D(farTL, farTR, SKYBLUE);
    DrawLine3D(farTR, farBR, SKYBLUE);
    DrawLine3D(farBR, farBL, SKYBLUE);
    DrawLine3D(farBL, farTL, SKYBLUE);
    DrawLine3D(nearTL, farTL, SKYBLUE);
    DrawLine3D(nearTR, farTR, SKYBLUE);
    DrawLine3D(nearBR, farBR, SKYBLUE);
    DrawLine3D(nearBL, farBL, SKYBLUE);
}

static void drawNearPlane(Camera3D *main, float aspect, float near, Mesh *mesh, Texture2D texture, float rotation)
{
    if (TEXTURE_MODE())
    {
        if (PERSPECTIVE_CORRECT())
        {
            buildNearPlaneQuad();
            updateNearPlaneQuad(main, aspect, near);
            nearPlaneQuadModel.materials[0].maps[MATERIAL_MAP_ALBEDO].texture = perspectiveCorrectTexture;
            rlEnableColorBlend();
            DrawModel(nearPlaneQuadModel, (Vector3){0, 0, 0}, 1.0f, WHITE);
        }
        else
        {
            perspectiveIncorrectTextureCapture(main, aspect, near, mesh, texture, rotation);
        }
    }
    else if (COLOR_MODE())
    {
        //TODO: reduce this please..
        flatColor(main, aspect, near, mesh, rotation);
    }
}

static void drawModelFilled(Model *model, Texture2D texture, float rotation)
{
    if (!(COLOR_MODE() || TEXTURE_MODE()))
        return;

    unsigned char *colorsBackup = model->meshes[0].colors;
    if (TEXTURE_MODE() && !COLOR_MODE())
        model->meshes[0].colors = NULL;

    model->materials[0].maps[MATERIAL_MAP_ALBEDO].texture.id = TEXTURE_MODE() ? texture.id : 0;
    DrawModelEx(*model, MODEL_POS, (Vector3){0, 1, 0}, RAD2DEG * rotation, MODEL_SCALE, WHITE);
    model->materials[0].maps[MATERIAL_MAP_ALBEDO].texture.id = 0;
    model->meshes[0].colors = colorsBackup;
}

static void drawWiresAndPoints(Model *model, float rotation)
{
    unsigned char *cacheColors = model->meshes[0].colors;
    model->meshes[0].colors = NULL;

    rlSetLineWidth(2.0f);
    DrawModelWiresEx(*model, MODEL_POS, (Vector3){0, 1, 0}, RAD2DEG * rotation, MODEL_SCALE, BLUE);
    rlSetPointSize(8.0f);
    DrawModelPointsEx(*model, MODEL_POS, (Vector3){0, 1, 0}, RAD2DEG * rotation, MODEL_SCALE, MAGENTA);

    model->meshes[0].colors = cacheColors;
}

static Vector3 nearPlaneIntersection(Camera3D *main, float near, Vector3 worldCoord)
{
    Vector3 viewDir = Vector3Normalize(Vector3Subtract(main->target, main->position));
    Vector3 mainCameraToPoint = Vector3Subtract(worldCoord, main->position);
    float depthAlongView = Vector3DotProduct(mainCameraToPoint, viewDir);
    if (depthAlongView <= 0.0f)
        return Vector3Add(main->position, Vector3Scale(viewDir, near));

    float scaleToNear = near / depthAlongView;
    return Vector3Add(main->position, Vector3Scale(mainCameraToPoint, scaleToNear));
}

static Vector3 translateRotateScale(Vector3 coordinate, Vector3 position, Vector3 scale, float rotation)
{
    Matrix matrix =
        MatrixMultiply(MatrixMultiply(MatrixScale(scale.x, scale.y, scale.z), MatrixRotateY(rotation)), MatrixTranslate(position.x, position.y, position.z));
    return Vector3Transform(coordinate, matrix);
}

static Vector3 inverseTranslateRotateScale(Vector3 coordinate, Vector3 position, Vector3 scale, float rotation)
{
    Matrix matrix =
        MatrixMultiply(MatrixMultiply(MatrixScale(scale.x, scale.y, scale.z), MatrixRotateY(rotation)), MatrixTranslate(position.x, position.y, position.z));
    Matrix inverseMatrix = MatrixInvert(matrix);
    return Vector3Transform(coordinate, inverseMatrix);
}

static Vector3 reflectPlane(Vector3 intersectionCoord, Camera3D *main, float near)
{
    Vector3 los, right, up;
    basisVector(main, &los, &right, &up);
    Vector3 centerNear = Vector3Add(main->position, Vector3Scale(los, near));
    Vector3 toClipPlane = Vector3Subtract(intersectionCoord, centerNear);
    float x = Vector3DotProduct(toClipPlane, right);
    float y = Vector3DotProduct(toClipPlane, up);
    return Vector3Add(centerNear, Vector3Add(Vector3Scale(right, x), Vector3Scale(up, -y)));
}

static void drawNearPlanePoints(Camera3D *main, float near, Model *nearPlanePointsModel, Mesh *mesh, float rotation, bool flipY)
{
    Mesh *nearPlanePointsMesh = &nearPlanePointsModel->meshes[0];
    int capacity = nearPlanePointsCapacity;
    int writtenVertices = 0;
    Vector3 los = Vector3Normalize(Vector3Subtract(main->target, main->position));

    for (int i = 0; i < mesh->triangleCount; i++)
    {
        int ia = mesh->indices ? mesh->indices[3 * i + 0] : 3 * i + 0;
        int ib = mesh->indices ? mesh->indices[3 * i + 1] : 3 * i + 1;
        int ic = mesh->indices ? mesh->indices[3 * i + 2] : 3 * i + 2;
        Vector3 a = (Vector3){mesh->vertices[3 * ia + 0], mesh->vertices[3 * ia + 1], mesh->vertices[3 * ia + 2]};
        Vector3 b = (Vector3){mesh->vertices[3 * ib + 0], mesh->vertices[3 * ib + 1], mesh->vertices[3 * ib + 2]};
        Vector3 c = (Vector3){mesh->vertices[3 * ic + 0], mesh->vertices[3 * ic + 1], mesh->vertices[3 * ic + 2]};
        a = translateRotateScale(a, MODEL_POS, MODEL_SCALE, rotation);
        b = translateRotateScale(b, MODEL_POS, MODEL_SCALE, rotation);
        c = translateRotateScale(c, MODEL_POS, MODEL_SCALE, rotation);

        //test if front facing or not
        Vector3 normal = Vector3Normalize(Vector3CrossProduct(Vector3Subtract(b, a), Vector3Subtract(c, a)));
        if (Vector3DotProduct(normal, los) > 0.0f)
            continue;

        Vector3 hits[3] = {nearPlaneIntersection(main, near, a), nearPlaneIntersection(main, near, b), nearPlaneIntersection(main, near, c)};

        for (int v = 0; v < 3 && writtenVertices < capacity; ++v)
        {
            Vector3 worldV = (Vector3[]){a, b, c}[v];
            Vector3 hit = flipY ? reflectPlane(hits[v], main, near) : hits[v];
            rlSetLineWidth(1.0f);
            DrawLine3D(worldV, hit, (Color){255, 0, 0, 80});
            nearPlanePointsMesh->vertices[3 * writtenVertices + 0] = hit.x;
            nearPlanePointsMesh->vertices[3 * writtenVertices + 1] = hit.y;
            nearPlanePointsMesh->vertices[3 * writtenVertices + 2] = hit.z;
            writtenVertices++;
        }
    }

    nearPlanePointsMesh->vertexCount = writtenVertices;

    rlSetPointSize(6.0f);
    DrawModelPoints(*nearPlanePointsModel, (Vector3){0, 0, 0}, 1.0f, MAGENTA);
}

//TODO: this has to exist somewhere in raylib math right???
static Vector3 rotatePointAboutAxis(Vector3 point, Vector3 axisPointA, Vector3 axisPointB, float angle)
{
    Vector3 axisDir = Vector3Normalize(Vector3Subtract(axisPointB, axisPointA));
    Vector3 localFromAxis = Vector3Subtract(point, axisPointA);
    Vector3 rotatedLocal = Vector3RotateByAxisAngle(localFromAxis, axisDir, angle);
    return Vector3Add(axisPointA, rotatedLocal);
}

static Vector3 worldToNDCCoord(float aspect, Camera3D *main, float near, float far, Vector3 worldCoord)
{
    Vector3 los, right, up;
    basisVector(main, &los, &right, &up);

    float halfFovy = DEG2RAD * main->fovy * 0.5f;
    float halfHNear = near * tanf(halfFovy);
    float halfWNear = halfHNear * aspect;

    float signedDepth = Vector3DotProduct(Vector3Subtract(worldCoord, main->position), los);
    if (signedDepth <= 0.0f)
    {
        return (Vector3){0.0f, 0.0f, 1.0f};
    }

    Vector3 intersectionCoord = nearPlaneIntersection(main, near, worldCoord);
    Vector3 centerNear = Vector3Add(main->position, Vector3Scale(los, near));
    Vector3 clipPlaneVector = Vector3Subtract(intersectionCoord, centerNear);

    float xNDC = Vector3DotProduct(clipPlaneVector, right) / halfWNear;
    float yNDC = Vector3DotProduct(clipPlaneVector, up) / halfHNear;

    float zNDC = ((far + near) - (2.0f * far * near) / signedDepth) / (far - near);

    return (Vector3){xNDC, yNDC, zNDC};
}

static Vector3 scaleNDCByNearClipPlane(Camera3D *main, Vector3 centerNear, float halfWNear, float halfHNear, float halfDepthNDCCube, Vector3 ndcCoord)
{
    Vector3 los, right, up;
    basisVector(main, &los, &right, &up);
    return Vector3Add(
        centerNear,
        Vector3Add(
            Vector3Scale(right, ndcCoord.x * halfWNear),
            Vector3Add(Vector3Scale(up, ndcCoord.y * halfHNear), Vector3Scale(los, ndcCoord.z * halfDepthNDCCube))));
}

static void worldToNDCSpace(Camera3D *main, float aspect, float near, float far, Model *world, Model *ndc, float rotation)
{
    Vector3 los = Vector3Normalize(Vector3Subtract(main->target, main->position));
    float halfFovy = DEG2RAD * main->fovy * 0.5f;
    float halfHNear = near * tanf(halfFovy);
    float halfWNear = ANISOTROPIC() ? (halfHNear * aspect) : halfHNear;
    float halfDepthNDC = ANISOTROPIC() ? (0.5f * (far - near)) : halfHNear;
    Vector3 centerNear = Vector3Add(main->position, Vector3Scale(los, near));
    Vector3 ndcCubeCenter = Vector3Add(centerNear, Vector3Scale(los, halfDepthNDC));
    drawNDCCube(main, ndcCubeCenter, halfWNear, halfHNear, halfDepthNDC);
    worldToNDCMesh(main, aspect, near, far, &world->meshes[0], &ndc->meshes[0], rotation);
}

static void worldToNDCMesh(Camera3D *main, float aspect, float near, float far, Mesh *world, Mesh *ndc, float rotation)
{
    Vector3 los = Vector3Normalize(Vector3Subtract(main->target, main->position));

    float halfFovy = DEG2RAD * main->fovy * 0.5f;
    float halfHNear = near * tanf(halfFovy);

    float halfWNear = ANISOTROPIC() ? (halfHNear * aspect) : halfHNear;
    float halfDepthNDC = ANISOTROPIC() ? (0.5f * (far - near)) : halfHNear;

    Vector3 centerNear = Vector3Add(main->position, Vector3Scale(los, near));
    Vector3 ndcCubeCenter = Vector3Add(centerNear, Vector3Scale(los, halfDepthNDC));
    float aspectForMapping = ANISOTROPIC() ? aspect : 1.0f;
    if (world->indices && ndc->indices)
    {
        memcpy(ndc->indices, world->indices, sizeof(unsigned short) * 3 * world->triangleCount);
        ndc->triangleCount = world->triangleCount;
    }

    for (int vertexIndex = 0; vertexIndex < world->vertexCount; vertexIndex++)
    {
        Vector3 objectVertex = {world->vertices[3 * vertexIndex + 0], world->vertices[3 * vertexIndex + 1], world->vertices[3 * vertexIndex + 2]};
        Vector3 worldVertex = translateRotateScale(objectVertex, MODEL_POS, MODEL_SCALE, rotation);
        Vector3 ndcCoord = worldToNDCCoord(aspectForMapping, main, near, far, worldVertex);
        Vector3 scaledNDCCoord = scaleNDCByNearClipPlane(main, ndcCubeCenter, halfWNear, halfHNear, halfDepthNDC, ndcCoord);
        Vector3 mappedObjectCoord = inverseTranslateRotateScale(scaledNDCCoord, MODEL_POS, MODEL_SCALE, rotation);

        ndc->vertices[3 * vertexIndex + 0] = mappedObjectCoord.x;
        ndc->vertices[3 * vertexIndex + 1] = mappedObjectCoord.y;
        ndc->vertices[3 * vertexIndex + 2] = mappedObjectCoord.z;
    }
}

static void buildNearPlaneQuad(void)
{
    //TODO: why is this having to be called every frame anyways? just make the quad and thats it i thought.
    if (nearPlaneQuad.vertexCount > 0 && nearPlaneQuadModel.meshCount > 0)
        return;

    nearPlaneQuad = (Mesh){0};
    nearPlaneQuad.vertexCount = 4;
    nearPlaneQuad.triangleCount = 2;
    nearPlaneQuad.vertices = (float *)RL_CALLOC(4 * 3, sizeof(float));
    nearPlaneQuad.texcoords = (float *)RL_CALLOC(4 * 2, sizeof(float));
    nearPlaneQuad.indices = (unsigned short *)RL_CALLOC(6, sizeof(unsigned short));
    float texcoords[8] = {0.0f, 0.0f, 1.0f, 0.0f, 1.0f, 1.0f, 0.0f, 1.0f};
    memcpy(nearPlaneQuad.texcoords, texcoords, sizeof(texcoords));
    unsigned short indices[6] = {0, 2, 1, 0, 3, 2};
    memcpy(nearPlaneQuad.indices, indices, sizeof(indices));
    nearPlaneQuadModel = LoadModelFromMesh(nearPlaneQuad);
}

static void updateNearPlaneQuad(Camera3D *main, float aspect, float near)
{
    //TODO: why is this having to be called every frame anyways? for aspect? why do we have to do this redundantly? just update only when needed?
    Vector3 los, right, up;
    basisVector(main, &los, &right, &up);
    float halfFovy = DEG2RAD * main->fovy * 0.5f;
    float halfH = near * tanf(halfFovy);
    float halfW = ANISOTROPIC() ? (halfH * aspect) : halfH;
    Vector3 centerNear = Vector3Add(main->position, Vector3Scale(los, near));
    Vector3 a = Vector3Add(centerNear, Vector3Add(Vector3Scale(up, +halfH), Vector3Scale(right, -halfW)));
    Vector3 b = Vector3Add(centerNear, Vector3Add(Vector3Scale(up, +halfH), Vector3Scale(right, +halfW)));
    Vector3 c = Vector3Add(centerNear, Vector3Add(Vector3Scale(up, -halfH), Vector3Scale(right, +halfW)));
    Vector3 d = Vector3Add(centerNear, Vector3Add(Vector3Scale(up, -halfH), Vector3Scale(right, -halfW)));

    float *v = nearPlaneQuad.vertices;
    v[0] = a.x;
    v[1] = a.y;
    v[2] = a.z;
    v[3] = b.x;
    v[4] = b.y;
    v[5] = b.z;
    v[6] = c.x;
    v[7] = c.y;
    v[8] = c.z;
    v[9] = d.x;
    v[10] = d.y;
    v[11] = d.z;
}

static Vector3 aspectCorrectNearPlane(Vector3 hit, Vector3 center, Vector3 right, Vector3 up, float xScale)
{
    //TODO: this is gross, and naming of functions is not clear and this is just bad and not clear at all in the chornology of the pipeline didactic
    Vector3 d = Vector3Subtract(hit, center);
    float xRight = Vector3DotProduct(d, right);
    float yUp = Vector3DotProduct(d, up);
    return Vector3Add(center, Vector3Add(Vector3Scale(right, xRight * xScale), Vector3Scale(up, yUp)));
}

static void perspectiveIncorrectTextureCapture(Camera3D *main, float aspect, float near, Mesh *mesh, Texture2D texture, float rotation)
{
    Vector3 los, right, up;
    basisVector(main, &los, &right, &up);
    Vector3 center = Vector3Add(main->position, Vector3Scale(los, near));
    float xScale = ANISOTROPIC() ? 1.0f : (1.0f / aspect);
    rlColor4f(1, 1, 1, 1); //TODO STILLL?????
    rlEnableTexture(texture.id);
    rlBegin(RL_TRIANGLES);

    for (int i = 0; i < mesh->triangleCount; i++)
    {
        int ia = mesh->indices ? mesh->indices[3 * i + 0] : 3 * i + 0;
        int ib = mesh->indices ? mesh->indices[3 * i + 1] : 3 * i + 1;
        int ic = mesh->indices ? mesh->indices[3 * i + 2] : 3 * i + 2;

        Vector3 a = (Vector3){mesh->vertices[3 * ia + 0], mesh->vertices[3 * ia + 1], mesh->vertices[3 * ia + 2]};
        Vector3 b = (Vector3){mesh->vertices[3 * ib + 0], mesh->vertices[3 * ib + 1], mesh->vertices[3 * ib + 2]};
        Vector3 c = (Vector3){mesh->vertices[3 * ic + 0], mesh->vertices[3 * ic + 1], mesh->vertices[3 * ic + 2]};

        a = translateRotateScale(a, MODEL_POS, MODEL_SCALE, rotation);
        b = translateRotateScale(b, MODEL_POS, MODEL_SCALE, rotation);
        c = translateRotateScale(c, MODEL_POS, MODEL_SCALE, rotation);

        Vector3 hita = nearPlaneIntersection(main, near, a);
        Vector3 hitb = nearPlaneIntersection(main, near, b);
        Vector3 hitc = nearPlaneIntersection(main, near, c);

        Vector3 pa = aspectCorrectNearPlane(hita, center, right, up, xScale);
        Vector3 pb = aspectCorrectNearPlane(hitb, center, right, up, xScale);
        Vector3 pc = aspectCorrectNearPlane(hitc, center, right, up, xScale);

        float sa = mesh->texcoords ? mesh->texcoords[2 * ia + 0] : 0.0f;
        float ta = mesh->texcoords ? mesh->texcoords[2 * ia + 1] : 0.0f;
        float sb = mesh->texcoords ? mesh->texcoords[2 * ib + 0] : 0.0f;
        float tb = mesh->texcoords ? mesh->texcoords[2 * ib + 1] : 0.0f;
        float sc = mesh->texcoords ? mesh->texcoords[2 * ic + 0] : 0.0f;
        float tc = mesh->texcoords ? mesh->texcoords[2 * ic + 1] : 0.0f;

        const unsigned char *ca = mesh->colors ? &mesh->colors[4 * ia] : NULL;
        const unsigned char *cb = mesh->colors ? &mesh->colors[4 * ib] : NULL;
        const unsigned char *cc = mesh->colors ? &mesh->colors[4 * ic] : NULL;
        //TODO: ewwwwwwwww
        if (COLOR_MODE() && ca)
            rlColor4ub(ca[0], ca[1], ca[2], ca[3]);
        else
            rlColor4ub(255, 255, 255, 255);
        rlTexCoord2f(sa, ta);
        rlVertex3f(pa.x, pa.y, pa.z);

        if (COLOR_MODE() && cb)
            rlColor4ub(cb[0], cb[1], cb[2], cb[3]);
        else
            rlColor4ub(255, 255, 255, 255);

        rlTexCoord2f(sb, tb);
        rlVertex3f(pb.x, pb.y, pb.z);

        if (COLOR_MODE() && cc)
            rlColor4ub(cc[0], cc[1], cc[2], cc[3]);
        else
            rlColor4ub(255, 255, 255, 255);
        rlTexCoord2f(sc, tc);
        rlVertex3f(pc.x, pc.y, pc.z);
    }

    rlEnd();
    rlDisableTexture();
}

//TODO: this LITERALLY does not work. color UVs and TEXCOORDS are NOT THE SAME (this is a good demonstration of why that is so i think????
static void flatColor(Camera3D *main, float aspect, float near, Mesh *mesh, float rotation)
{
    Vector3 los, right, up;
    basisVector(main, &los, &right, &up);
    Vector3 centerNear = Vector3Add(main->position, Vector3Scale(los, near));
    float xScale = ANISOTROPIC() ? 1.0f : (1.0f / aspect);
    rlDisableTexture();
    rlBegin(RL_TRIANGLES);

    for (int i = 0; i < mesh->triangleCount; i++)
    {
        int ia = mesh->indices ? mesh->indices[3 * i + 0] : 3 * i + 0;
        int ib = mesh->indices ? mesh->indices[3 * i + 1] : 3 * i + 1;
        int ic = mesh->indices ? mesh->indices[3 * i + 2] : 3 * i + 2;

        Vector3 a = (Vector3){mesh->vertices[3 * ia + 0], mesh->vertices[3 * ia + 1], mesh->vertices[3 * ia + 2]};
        Vector3 b = (Vector3){mesh->vertices[3 * ib + 0], mesh->vertices[3 * ib + 1], mesh->vertices[3 * ib + 2]};
        Vector3 c = (Vector3){mesh->vertices[3 * ic + 0], mesh->vertices[3 * ic + 1], mesh->vertices[3 * ic + 2]};

        a = translateRotateScale(a, MODEL_POS, MODEL_SCALE, rotation);
        b = translateRotateScale(b, MODEL_POS, MODEL_SCALE, rotation);
        c = translateRotateScale(c, MODEL_POS, MODEL_SCALE, rotation);

        Vector3 hita = nearPlaneIntersection(main, near, a);
        Vector3 hitb = nearPlaneIntersection(main, near, b);
        Vector3 hitc = nearPlaneIntersection(main, near, c);

        Vector3 pa = aspectCorrectNearPlane(hita, centerNear, right, up, xScale);
        Vector3 pb = aspectCorrectNearPlane(hitb, centerNear, right, up, xScale);
        Vector3 pc = aspectCorrectNearPlane(hitc, centerNear, right, up, xScale);

        const unsigned char *ca = &mesh->colors[4 * ia];
        const unsigned char *cb = &mesh->colors[4 * ib];
        const unsigned char *cc = &mesh->colors[4 * ic];

        rlColor4ub(ca[0], ca[1], ca[2], ca[3]);
        rlVertex3f(pa.x, pa.y, pa.z);
        rlColor4ub(cb[0], cb[1], cb[2], cb[3]);
        rlVertex3f(pb.x, pb.y, pb.z);
        rlColor4ub(cc[0], cc[1], cc[2], cc[3]);
        rlVertex3f(pc.x, pc.y, pc.z);
    }

    rlEnd();
    rlEnableTexture(rlGetTextureIdDefault());
}

static void applyBarycentricPalette(Mesh *mesh)
{
    if (mesh->colors == NULL)
        //TODO: oof
        mesh->colors = (unsigned char *)RL_CALLOC(mesh->vertexCount * 4, sizeof(unsigned char));
    for (int i = 0; i < mesh->triangleCount; i++)
    {
        int ia = mesh->indices ? mesh->indices[3 * i + 0] : 3 * i + 0;
        int ib = mesh->indices ? mesh->indices[3 * i + 1] : 3 * i + 1;
        int ic = mesh->indices ? mesh->indices[3 * i + 2] : 3 * i + 2;

        //TODO: can we more cleanly achieve this? like i just want to set a color like = RED = GREEN = BLUE, i dont want ot fuck with pointer offsets and arrays...
        mesh->colors[4 * ia + 0] = 255;
        mesh->colors[4 * ia + 1] = 0;
        mesh->colors[4 * ia + 2] = 0;
        mesh->colors[4 * ia + 3] = 255;

        mesh->colors[4 * ib + 0] = 0;
        mesh->colors[4 * ib + 1] = 255;
        mesh->colors[4 * ib + 2] = 0;
        mesh->colors[4 * ib + 3] = 255;

        mesh->colors[4 * ic + 0] = 0;
        mesh->colors[4 * ic + 1] = 0;
        mesh->colors[4 * ic + 2] = 255;
        mesh->colors[4 * ic + 3] = 255;
    }
}

static void perspectiveCorrectTextureCapture(Camera3D *main, Model *model, Texture2D texture, float rotation)
{
    unsigned char *cacheVertexColors = model->meshes[0].colors;
    if (TEXTURE_MODE() && !COLOR_MODE())
        model->meshes[0].colors = NULL;

    //TODO: still a bit confused about this one for maintaining the cohesion between the near plane and the actual model visual
    ClearBackground(BLACK);
    BeginMode3D(*main);
    Texture2D previousTexture = model->materials[0].maps[MATERIAL_MAP_ALBEDO].texture;
    model->materials[0].maps[MATERIAL_MAP_ALBEDO].texture = texture;
    DrawModelEx(*model, MODEL_POS, (Vector3){0, 1, 0}, RAD2DEG * rotation, MODEL_SCALE, WHITE);
    model->materials[0].maps[MATERIAL_MAP_ALBEDO].texture = previousTexture;
    EndMode3D();

    Image rgba = LoadImageFromScreen();
    ImageFormat(&rgba, PIXELFORMAT_UNCOMPRESSED_R8G8B8A8);
    if (TEXTURE_MODE() && !COLOR_MODE())
        model->meshes[0].colors = cacheVertexColors;

    ClearBackground(BLACK);
    BeginMode3D(*main);
    Texture2D cacheTexture = model->materials[0].maps[MATERIAL_MAP_ALBEDO].texture;
    Color cacheMaterialColor = model->materials[0].maps[MATERIAL_MAP_ALBEDO].color;
    model->materials[0].maps[MATERIAL_MAP_ALBEDO].texture = (Texture2D){0};
    model->materials[0].maps[MATERIAL_MAP_ALBEDO].color = WHITE;
    DrawModelEx(*model, MODEL_POS, (Vector3){0, 1, 0}, RAD2DEG * rotation, MODEL_SCALE, WHITE);
    model->materials[0].maps[MATERIAL_MAP_ALBEDO].texture = cacheTexture;
    model->materials[0].maps[MATERIAL_MAP_ALBEDO].color = cacheMaterialColor;
    EndMode3D();

    Image mask = LoadImageFromScreen();
    punchOutAlphaFromImage(&rgba, &mask, 1);

    if (perspectiveCorrectTexture.id != 0) //TODO: the fuck ew??
        UpdateTexture(perspectiveCorrectTexture, rgba.data);
    else
        perspectiveCorrectTexture = LoadTextureFromImage(rgba);

    UnloadImage(mask);
    UnloadImage(rgba);
}

static void punchOutAlphaFromImage(Image *rgba, const Image *mask, unsigned char threshold)
{
    Image maskCopy = ImageCopy(*mask);
    ImageFormat(&maskCopy, PIXELFORMAT_UNCOMPRESSED_R8G8B8A8);
    ImageFormat(rgba, PIXELFORMAT_UNCOMPRESSED_R8G8B8A8);
    unsigned char *rgbaPixels = rgba->data;
    unsigned char *maskPixels = maskCopy.data;
    int pixelCount = rgba->width * rgba->height;
    for (int i = 0; i < pixelCount; ++i)
    {
        unsigned int red = maskPixels[4 * i + 0];
        unsigned int green = maskPixels[4 * i + 1];
        unsigned int blue = maskPixels[4 * i + 2];
        unsigned int luma = (red * 30u + green * 59u + blue * 11u) / 100u;
        rgbaPixels[4 * i + 3] = (luma > threshold) ? 255 : 0;
    }
    UnloadImage(maskCopy);
}

static void moveJugemuOrbital(Camera3D *jugemu, float deltaTime)
{
    float radius = Vector3Length(jugemu->position);
    float azimuth = atan2f(jugemu->position.z, jugemu->position.x);
    float horizontalRadius = sqrtf(jugemu->position.x * jugemu->position.x + jugemu->position.z * jugemu->position.z);
    float elevation = atan2f(jugemu->position.y, horizontalRadius);
    const float LONG_SPEED = 1.5f;
    const float LAT_SPEED = 1.0f;
    const float ZOOM_SPEED = 2.0f;
    const float ROLL_SPEED = 2.0f;
    float rollDeltaRadians = 0.0f;

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
    elevation = Clamp(elevation, -PI * 0.5f + 0.0001f, PI * 0.5f - 0.0001f);
    jugemu->position.x = radius * cosf(elevation) * cosf(azimuth);
    jugemu->position.y = radius * sinf(elevation);
    jugemu->position.z = radius * cosf(elevation) * sinf(azimuth);
    Vector3 los = Vector3Normalize(Vector3Subtract((Vector3){0, 0, 0}, jugemu->position));
    Vector3 up = rotatePointAboutAxis(jugemu->up, (Vector3){0, 0, 0}, los, rollDeltaRadians);
    jugemu->target = (Vector3){0, 0, 0};
    jugemu->up = Vector3Normalize(up);
}