#include "raylib.h"
#include "raymath.h"
#include "rlgl.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

static const int WIDTH = 640 * 1.5;
static const int HEIGHT = 480 * 1.5;
static const float ANGULAR_VELOCITY = 1.25f;
static const float FOVY_PERSPECTIVE = 60.0f;
static const Vector3 MODEL_POS = {0.0f, 0.0f, 0.0f};
static const Vector3 MODEL_SCALE = {1.0f, 1.0f, 1.0f};
static const Vector3 MAIN_POS = {0.0f, 0.0f, 2.0f};
static const Vector3 JUGEMU_POS_ISO = {3.0f, 1.0f, 3.0f};

enum
{
    FLAG_NDC = 1u << 0,                 // N
    FLAG_ASPECT = 1u << 1,              // Q
    FLAG_PERSPECTIVE_CORRECT = 1u << 2, // P
    FLAG_PAUSE = 1u << 3,               // SPACE
    FLAG_COLOR_MODE = 1u << 4,          // C
    FLAG_TEXTURE_MODE = 1u << 5         // T
};
static unsigned int gFlags = FLAG_ASPECT | FLAG_COLOR_MODE;

#define NDC_SPACE() ((gFlags & FLAG_NDC) != 0)
#define ANISOTROPIC() ((gFlags & FLAG_ASPECT) != 0)
#define PERSPECTIVE_CORRECT() ((gFlags & FLAG_PERSPECTIVE_CORRECT) != 0)
#define PAUSED() ((gFlags & FLAG_PAUSE) != 0)
#define COLOR_MODE() ((gFlags & FLAG_COLOR_MODE) != 0)
#define TEXTURE_MODE() ((gFlags & FLAG_TEXTURE_MODE) != 0)
#define TOGGLE(K, F)                                                                                                                                           \
    do                                                                                                                                                         \
    {                                                                                                                                                          \
        if (IsKeyPressed(K)) gFlags ^= (F);                                                                                                                    \
    } while (0)

static void drawBasisVector(Camera3D *main);
static void updateSpatialWire(Camera3D *main, float aspect, float near, float far, Mesh *spatialWire);
static void drawSpatialWire(Mesh spatialWire);
static void showOnlyFrontFace(Mesh *mesh);
static void drawModelFilled(Model *model, Texture2D texture, float rotation);
static void drawWiresAndPoints(Model *model, float rotation);

static int nearPlanePointsCapacity = 0; //TODO: does this need to be global really?
static void drawNearPlanePoints(Camera3D *main, float near, Model *nearPlanePointsModel, Mesh *mesh, float rotation, bool flipY);
static void drawNearPlane(
    Camera3D *main, float aspect, float near, Model spatialWireModel, Mesh *mesh, Texture2D meshTexture, Texture2D perspectiveCorrectTexture, float rotation);

static void worldToNDCSpace(Camera3D *main, float aspect, float near, float far, Model *world, Model *ndc, float rotation);

static void perspectiveIncorrectCapture(Camera3D *main, float aspect, float near, Mesh *mesh, Texture2D meshTexture, float rotation);
static void perspectiveCorrectCapture(Camera3D *main, Model *model, Texture2D meshTexture, Texture2D *perspectiveCorrectTexture, float rotation);
static void alphaMaskPunchOut(Image *rgba, const Image *mask, unsigned char threshold);
static void fillVertexColors(Mesh *mesh);
static void moveJugemuOrbital(Camera3D *jugemu, float deltaTime);

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
    Model worldModel = LoadModelFromMesh(GenMeshCube(1.0f, 1.0f, 1.0f));
    // Model worldModel = LoadModelFromMesh(GenMeshKnot(1.0f, 1.0f, 16, 128));

    Mesh ndcMesh = (Mesh){0};
    ndcMesh.vertexCount = worldModel.meshes[0].vertexCount;
    ndcMesh.triangleCount = worldModel.meshes[0].triangleCount;
    ndcMesh.vertices = RL_CALLOC(ndcMesh.vertexCount * 3, sizeof(float));
    ndcMesh.texcoords = RL_CALLOC(ndcMesh.vertexCount * 2, sizeof(float));
    ndcMesh.indices = RL_CALLOC(ndcMesh.triangleCount * 3, sizeof(unsigned short));

    if (worldModel.meshes[0].texcoords) memcpy(ndcMesh.texcoords, worldModel.meshes[0].texcoords, sizeof(float) * 2 * ndcMesh.vertexCount);
    if (worldModel.meshes[0].indices) memcpy(ndcMesh.indices, worldModel.meshes[0].indices, sizeof(unsigned short) * 3 * ndcMesh.triangleCount);

    Model ndcModel = LoadModelFromMesh(ndcMesh);

    Mesh nearPlanePoints = (Mesh){0};
    nearPlanePoints.vertexCount = 3 * worldModel.meshes[0].triangleCount;
    nearPlanePointsCapacity = nearPlanePoints.vertexCount;
    nearPlanePoints.vertices = RL_CALLOC(nearPlanePoints.vertexCount * 3, sizeof(float));
    Model nearPlanePointsModel = LoadModelFromMesh(nearPlanePoints);

    Image checkeredImage = GenImageChecked(16, 16, 4, 4, BLACK, WHITE);
    Texture2D checkeredTexture = LoadTextureFromImage(checkeredImage);
    UnloadImage(checkeredImage);

    worldModel.materials[0].maps[MATERIAL_MAP_ALBEDO].texture = checkeredTexture;
    fillVertexColors(&worldModel.meshes[0]);
    ndcModel.materials[0].maps[MATERIAL_MAP_ALBEDO].texture = checkeredTexture;
    fillVertexColors(&ndcModel.meshes[0]);

    Texture2D perspectiveCorrectTexture = {0}; //TODO: still wacky to do this and pass by reference everywhere i dont like it... its not neccessary
    Mesh spatialWire = GenMeshCube(1.0f, 1.0f, 1.0f);
    showOnlyFrontFace(&spatialWire);
    Model spatialWireModel = LoadModelFromMesh(spatialWire);
    while (!WindowShouldClose())
    {
        aspect = (float)GetScreenWidth() / (float)GetScreenHeight();
        TOGGLE(KEY_N, FLAG_NDC);
        TOGGLE(KEY_Q, FLAG_ASPECT);
        TOGGLE(KEY_P, FLAG_PERSPECTIVE_CORRECT);
        TOGGLE(KEY_SPACE, FLAG_PAUSE);
        TOGGLE(KEY_C, FLAG_COLOR_MODE);
        TOGGLE(KEY_T, FLAG_TEXTURE_MODE);

        moveJugemuOrbital(&jugemu, GetFrameTime());
        if (!PAUSED()) meshRotation -= ANGULAR_VELOCITY * GetFrameTime();
        Model *capturedModel = NDC_SPACE() ? &ndcModel : &worldModel;
        Mesh *capturedMesh = NDC_SPACE() ? &ndcModel.meshes[0] : &worldModel.meshes[0];
        if (NDC_SPACE()) worldToNDCSpace(&main, aspect, near, far, &worldModel, &ndcModel, meshRotation);

        BeginDrawing();
        if (PERSPECTIVE_CORRECT() && TEXTURE_MODE())
            perspectiveCorrectCapture(&main, capturedModel, checkeredTexture, &perspectiveCorrectTexture, meshRotation);
        ClearBackground(BLACK);

        BeginMode3D(jugemu);
        drawBasisVector(&main);
        updateSpatialWire(&main, ANISOTROPIC() ? aspect : 1.0f, near, far, &spatialWire);
        drawSpatialWire(spatialWire);
        drawModelFilled(capturedModel, checkeredTexture, meshRotation);
        drawWiresAndPoints(capturedModel, meshRotation);
        drawNearPlanePoints(&main, near, &nearPlanePointsModel, capturedMesh, meshRotation, NDC_SPACE());
        drawNearPlane(&main, ANISOTROPIC() ? aspect : 1.0f, near, spatialWireModel, capturedMesh, checkeredTexture, perspectiveCorrectTexture, meshRotation);
        EndMode3D();

        DrawText("Arrows/WASD: move camera | R: reset camera | Spacebar: pause", 12, 40, 20, WHITE);
        DrawText("E: NDC | Q: ASPECT | P: PERSPECTIVE | C: COLOR | T: TEXTURE", 12, 14, 20, WHITE);
        DrawText("SPACE:", 12, 66, 20, WHITE);
        DrawText(NDC_SPACE() ? "NDC" : "WORLD", 180, 66, 20, WHITE);
        DrawText("ASPECT:", 12, 92, 20, WHITE);
        DrawText(ANISOTROPIC() ? "ANISOTROPIC" : "ISOTROPIC", 180, 92, 20, WHITE);
        DrawText("TEXTURE:", 12, 118, 20, WHITE);
        DrawText(TEXTURE_MODE() ? "ON" : "OFF", 180, 118, 20, WHITE);
        DrawText("PERSPECTIVE:", 12, 144, 20, WHITE);
        DrawText(PERSPECTIVE_CORRECT() ? "CORRECT" : "INCORRECT", 180, 144, 20, WHITE);
        DrawText("COLORS:", 12, 170, 20, WHITE);
        DrawText(COLOR_MODE() ? "ON" : "OFF", 180, 170, 20, WHITE);
        EndDrawing();
    }

    UnloadModel(worldModel);
    UnloadModel(ndcModel);
    UnloadModel(nearPlanePointsModel);
    UnloadModel(spatialWireModel);
    UnloadTexture(checkeredTexture);
    CloseWindow();
    return 0;
}

static void showOnlyFrontFace(Mesh *mesh)
{
    if (!mesh->colors) mesh->colors = (unsigned char *)RL_CALLOC(mesh->vertexCount, sizeof(Color));
    Color *colors = (Color *)mesh->colors;
    for (int i = 0; i < mesh->vertexCount; i++)
        colors[i] = (Color){255, 255, 255, 0};
    //first 4 indices is the front face of the mesh
    for (int i = 0; i < 4; i++)
    {
        colors[i].a = 255;
    }
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

static void drawBasisVector(Camera3D *main)
{
    Vector3 los, right, up;
    basisVector(main, &los, &right, &up);
    DrawLine3D(main->position, Vector3Add(main->position, right), PURPLE);
    DrawLine3D(main->position, Vector3Add(main->position, up), DARKGREEN);
    DrawLine3D(main->position, Vector3Add(main->position, los), SKYBLUE);
}

static void updateSpatialWire(Camera3D *main, float aspect, float near, float far, Mesh *spatialWire)
{
    Vector3 los, right, up;
    basisVector(main, &los, &right, &up);
    float halfFovy = DEG2RAD * main->fovy * 0.5f;
    float halfHNear = near * tanf(halfFovy);
    float halfWNear = ANISOTROPIC() ? (halfHNear * aspect) : halfHNear;
    float halfHFar = far * tanf(halfFovy);
    float halfWFar = ANISOTROPIC() ? (halfHFar * aspect) : halfHFar;
    Vector3 centerNear = Vector3Add(main->position, Vector3Scale(los, near));
    float halfDepth = NDC_SPACE() ? (ANISOTROPIC() ? (0.5f * (far - near)) : halfHNear) : (0.5f * (far - near));
    float depthSpan = 2.0f * halfDepth;
    float farHalfW = NDC_SPACE() ? halfWNear : halfWFar;
    float farHalfH = NDC_SPACE() ? halfHNear : halfHFar;

    for (int i = 0; i < spatialWire->vertexCount; ++i)
    {
        Vector3 vertex = (Vector3){spatialWire->vertices[3 * i + 0], spatialWire->vertices[3 * i + 1], spatialWire->vertices[3 * i + 2]};
        Vector3 offset = Vector3Subtract(vertex, centerNear);
        float xSign = (Vector3DotProduct(offset, right) >= 0.0f) ? 1.0f : -1.0f;
        float ySign = (Vector3DotProduct(offset, up) >= 0.0f) ? 1.0f : -1.0f;
        float farMask = (Vector3DotProduct(offset, los) > halfDepth) ? 1.0f : 0.0f;
        float finalHalfWidth = halfWNear + farMask * (farHalfW - halfWNear);
        float finalHalfHeight = halfHNear + farMask * (farHalfH - halfHNear);
        Vector3 faceCenter = Vector3Add(centerNear, Vector3Scale(los, farMask * depthSpan));
        Vector3 result = Vector3Add(faceCenter, Vector3Add(Vector3Scale(right, xSign * finalHalfWidth), Vector3Scale(up, ySign * finalHalfHeight)));
        spatialWire->vertices[3 * i + 0] = result.x;
        spatialWire->vertices[3 * i + 1] = result.y;
        spatialWire->vertices[3 * i + 2] = result.z;
    }
}

static void drawSpatialWire(Mesh spatialWire)
{
    //DrawModelWiresEx(spatialWireModel, MODEL_POS, (Vector3){0, 1, 0}, 0.0f, MODEL_SCALE, WHITE);
    static const int frontFace[4][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
    static const int backFace[4][2] = {{4, 5}, {5, 6}, {6, 7}, {7, 4}};
    static const int connectingFaces[4][2] = {{0, 4}, {1, 7}, {2, 6}, {3, 5}};
    static const int (*faces[3])[2] = {frontFace, backFace, connectingFaces};

    rlSetLineWidth(2.0f);
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 4; j++)
        {
            const int (*edges)[2] = faces[i];
            int edgeIndexA = edges[j][0] * 3;
            int edgeIndexB = edges[j][1] * 3;
            Vector3 startPosition = {spatialWire.vertices[edgeIndexA + 0], spatialWire.vertices[edgeIndexA + 1], spatialWire.vertices[edgeIndexA + 2]};
            Vector3 endPosition = {spatialWire.vertices[edgeIndexB + 0], spatialWire.vertices[edgeIndexB + 1], spatialWire.vertices[edgeIndexB + 2]};
            DrawLine3D(startPosition, endPosition, WHITE);
        }
}

static void drawModelFilled(Model *model, Texture2D texture, float rotation)
{
    if (!(COLOR_MODE() || TEXTURE_MODE()))
    {
        return;
    }
    unsigned char *colorsBackup = model->meshes[0].colors;
    if (TEXTURE_MODE() && !COLOR_MODE()) model->meshes[0].colors = NULL;
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
    // rlSetPointSize(8.0f);
    // DrawModelPointsEx(*model, MODEL_POS, (Vector3){0, 1, 0}, RAD2DEG * rotation, MODEL_SCALE, MAGENTA);
    model->meshes[0].colors = cacheColors;
}

static Vector3 nearPlaneIntersection(Camera3D *main, float near, Vector3 worldCoord)
{
    Vector3 viewDir = Vector3Normalize(Vector3Subtract(main->target, main->position));
    Vector3 mainCameraToPoint = Vector3Subtract(worldCoord, main->position);
    float depthAlongView = Vector3DotProduct(mainCameraToPoint, viewDir);
    if (depthAlongView <= 0.0f)
    {
        return Vector3Add(main->position, Vector3Scale(viewDir, near));
    }
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
    Vector3 centerNearPlane = Vector3Add(main->position, Vector3Scale(los, near));
    Vector3 toClipPlane = Vector3Subtract(intersectionCoord, centerNearPlane);
    float x = Vector3DotProduct(toClipPlane, right);
    float y = Vector3DotProduct(toClipPlane, up);
    return Vector3Add(centerNearPlane, Vector3Add(Vector3Scale(right, x), Vector3Scale(up, -y)));
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
        if (Vector3DotProduct(Vector3Normalize(Vector3CrossProduct(Vector3Subtract(b, a), Vector3Subtract(c, a))), los) > 0.0f)
        {
            continue;
        }

        Vector3 hits[3] = {nearPlaneIntersection(main, near, a), nearPlaneIntersection(main, near, b), nearPlaneIntersection(main, near, c)};

        for (int j = 0; j < 3 && writtenVertices < capacity; ++j)
        {
            Vector3 worldV = (Vector3[]){a, b, c}[j];
            Vector3 hit = flipY ? reflectPlane(hits[j], main, near) : hits[j];
            rlSetLineWidth(1.0f);
            DrawLine3D(worldV, hit, (Color){255, 0, 0, 80});
            nearPlanePointsMesh->vertices[3 * writtenVertices + 0] = hit.x;
            nearPlanePointsMesh->vertices[3 * writtenVertices + 1] = hit.y;
            nearPlanePointsMesh->vertices[3 * writtenVertices + 2] = hit.z;
            writtenVertices++;
        }
    }

    nearPlanePointsMesh->vertexCount = writtenVertices;
    // rlSetPointSize(6.0f);
    // DrawModelPoints(*nearPlanePointsModel, MODEL_POS, 1.0f, MAGENTA);
}

static void drawNearPlane(
    Camera3D *main, float aspect, float near, Model spatialWireModel, Mesh *mesh, Texture2D meshTexture, Texture2D perspectiveCorrectTexture, float rotation)
{
    if (PERSPECTIVE_CORRECT() && TEXTURE_MODE())
    {
        spatialWireModel.materials[0].maps[MATERIAL_MAP_ALBEDO].texture = perspectiveCorrectTexture;
        DrawModel(spatialWireModel, MODEL_POS, 1.0f, WHITE);
        return;
    }
    perspectiveIncorrectCapture(main, aspect, near, mesh, meshTexture, rotation);
}

//TODO: this has to exist somewhere in raylib math right???
static Vector3 rotatePointAboutAxis(Vector3 point, Vector3 axisPointA, Vector3 axisPointB, float angle)
{
    Vector3 axisDir = Vector3Normalize(Vector3Subtract(axisPointB, axisPointA));
    Vector3 localFromAxis = Vector3Subtract(point, axisPointA);
    Vector3 rotatedLocal = Vector3RotateByAxisAngle(localFromAxis, axisDir, angle);
    return Vector3Add(axisPointA, rotatedLocal);
}

static void worldToNDCSpace(Camera3D *main, float aspect, float near, float far, Model *world, Model *ndc, float rotation)
{
    Vector3 los, right, up;
    basisVector(main, &los, &right, &up);
    float halfFovy = DEG2RAD * main->fovy * 0.5f;
    float halfHNear = near * tanf(halfFovy);
    float halfWNear = ANISOTROPIC() ? (halfHNear * aspect) : halfHNear;
    float halfDepthNDC = ANISOTROPIC() ? (0.5f * (far - near)) : halfHNear;
    Vector3 centerNearPlane = Vector3Add(main->position, Vector3Scale(los, near));
    Vector3 centerNDCCube = Vector3Add(centerNearPlane, Vector3Scale(los, halfDepthNDC));
    Mesh *worldMesh = &world->meshes[0];
    Mesh *ndcMesh = &ndc->meshes[0];

    for (int i = 0; i < worldMesh->vertexCount; i++)
    {
        Vector3 objectVertex = {worldMesh->vertices[3 * i + 0], worldMesh->vertices[3 * i + 1], worldMesh->vertices[3 * i + 2]};
        Vector3 worldVertex = translateRotateScale(objectVertex, MODEL_POS, MODEL_SCALE, rotation);
        float signedDepth = Vector3DotProduct(Vector3Subtract(worldVertex, main->position), los);
        Vector3 ndcCoord;
        if (signedDepth <= 0.0f)
        {
            ndcCoord = (Vector3){0.0f, 0.0f, 1.0f};
        }
        else
        {
            Vector3 intersectionCoord = nearPlaneIntersection(main, near, worldVertex);
            Vector3 clipPlaneVector = Vector3Subtract(intersectionCoord, centerNearPlane);
            float xNDC = Vector3DotProduct(clipPlaneVector, right) / halfWNear;
            float yNDC = Vector3DotProduct(clipPlaneVector, up) / halfHNear;
            float zNDC = ((far + near) - (2.0f * far * near) / signedDepth) / (far - near);
            ndcCoord = (Vector3){xNDC, yNDC, zNDC};
        }
        //TODO: ugly as hell, clean it up somehow PLEASE
        Vector3 scaledNDCCoord = Vector3Add(
            centerNDCCube,
            Vector3Add(
                Vector3Scale(right, ndcCoord.x * halfWNear),
                Vector3Add(Vector3Scale(up, ndcCoord.y * halfHNear), Vector3Scale(los, ndcCoord.z * halfDepthNDC))));
        Vector3 mappedObjectCoord = inverseTranslateRotateScale(scaledNDCCoord, MODEL_POS, MODEL_SCALE, rotation);
        ndcMesh->vertices[3 * i + 0] = mappedObjectCoord.x;
        ndcMesh->vertices[3 * i + 1] = mappedObjectCoord.y;
        ndcMesh->vertices[3 * i + 2] = mappedObjectCoord.z;
    }
}

static Vector3 aspectCorrectNearPlane(Vector3 hit, Vector3 center, Vector3 right, Vector3 up, float xScale)
{
    Vector3 centerDistance = Vector3Subtract(hit, center);
    float xRight = Vector3DotProduct(centerDistance, right);
    float yUp = Vector3DotProduct(centerDistance, up);
    return Vector3Add(center, Vector3Add(Vector3Scale(right, xRight * xScale), Vector3Scale(up, yUp)));
}

static void perspectiveIncorrectCapture(Camera3D *main, float aspect, float near, Mesh *mesh, Texture2D meshTexture, float rotation)
{
    Vector3 los, right, up;
    basisVector(main, &los, &right, &up);
    Vector3 centerNearPlane = Vector3Add(main->position, Vector3Scale(los, near));
    float xScale = ANISOTROPIC() ? 1.0f : 1.0f / aspect;
    rlColor4f(WHITE.r, WHITE.g, WHITE.b, WHITE.a);
    if (TEXTURE_MODE())
        rlEnableTexture(meshTexture.id);
    else
        rlDisableTexture();

    if (!TEXTURE_MODE() && !COLOR_MODE())
    {
        rlEnableWireMode();
        rlColor4ub(BLUE.r, BLUE.g, BLUE.b, BLUE.a);
    }
    rlBegin(RL_TRIANGLES);

    for (int i = 0; i < mesh->triangleCount; i++)
    {
        int ia = mesh->indices ? mesh->indices[3 * i + 0] : 3 * i + 0;
        int ib = mesh->indices ? mesh->indices[3 * i + 1] : 3 * i + 1;
        int ic = mesh->indices ? mesh->indices[3 * i + 2] : 3 * i + 2;
        //TODO: the vertex mutation stages is gross here idk how to clean it up though
        Vector3 vertexA = (Vector3){mesh->vertices[3 * ia + 0], mesh->vertices[3 * ia + 1], mesh->vertices[3 * ia + 2]};
        Vector3 vertexB = (Vector3){mesh->vertices[3 * ib + 0], mesh->vertices[3 * ib + 1], mesh->vertices[3 * ib + 2]};
        Vector3 vertexC = (Vector3){mesh->vertices[3 * ic + 0], mesh->vertices[3 * ic + 1], mesh->vertices[3 * ic + 2]};

        vertexA = translateRotateScale(vertexA, MODEL_POS, MODEL_SCALE, rotation);
        vertexB = translateRotateScale(vertexB, MODEL_POS, MODEL_SCALE, rotation);
        vertexC = translateRotateScale(vertexC, MODEL_POS, MODEL_SCALE, rotation);

        vertexA = aspectCorrectNearPlane(nearPlaneIntersection(main, near, vertexA), centerNearPlane, right, up, xScale);
        vertexB = aspectCorrectNearPlane(nearPlaneIntersection(main, near, vertexB), centerNearPlane, right, up, xScale);
        vertexC = aspectCorrectNearPlane(nearPlaneIntersection(main, near, vertexC), centerNearPlane, right, up, xScale);

        if (COLOR_MODE()) rlColor4ub(mesh->colors[4 * ia + 0], mesh->colors[4 * ia + 1], mesh->colors[4 * ia + 2], mesh->colors[4 * ia + 3]);
        if (TEXTURE_MODE()) rlTexCoord2f(mesh->texcoords[2 * ia + 0], mesh->texcoords[2 * ia + 1]);
        rlVertex3f(vertexA.x, vertexA.y, vertexA.z);

        if (COLOR_MODE()) rlColor4ub(mesh->colors[4 * ib + 0], mesh->colors[4 * ib + 1], mesh->colors[4 * ib + 2], mesh->colors[4 * ib + 3]);
        if (TEXTURE_MODE()) rlTexCoord2f(mesh->texcoords[2 * ib + 0], mesh->texcoords[2 * ib + 1]);
        rlVertex3f(vertexB.x, vertexB.y, vertexB.z);

        if (COLOR_MODE()) rlColor4ub(mesh->colors[4 * ic + 0], mesh->colors[4 * ic + 1], mesh->colors[4 * ic + 2], mesh->colors[4 * ic + 3]);
        if (TEXTURE_MODE()) rlTexCoord2f(mesh->texcoords[2 * ic + 0], mesh->texcoords[2 * ic + 1]);
        rlVertex3f(vertexC.x, vertexC.y, vertexC.z);
    }

    rlEnd();
    rlDisableTexture();
    rlDisableWireMode();
}

static void perspectiveCorrectCapture(Camera3D *main, Model *model, Texture2D meshTexture, Texture2D *perspectiveCorrectTexture, float rotation)
{
    unsigned char *cacheVertexColors = model->meshes[0].colors;
    if (TEXTURE_MODE() && !COLOR_MODE()) model->meshes[0].colors = NULL;
    ClearBackground(BLACK);
    BeginMode3D(*main);
    //TODO: this is still needded to avoid the model in space and model proejected to near plane texture divergence... ugly though
    Texture2D previousTexture = model->materials[0].maps[MATERIAL_MAP_ALBEDO].texture;
    model->materials[0].maps[MATERIAL_MAP_ALBEDO].texture = meshTexture;
    DrawModelEx(*model, MODEL_POS, (Vector3){0, 1, 0}, RAD2DEG * rotation, MODEL_SCALE, WHITE);
    model->materials[0].maps[MATERIAL_MAP_ALBEDO].texture = previousTexture;
    EndMode3D();
    Image rgba = LoadImageFromScreen();
    ImageFormat(&rgba, PIXELFORMAT_UNCOMPRESSED_R8G8B8A8);
    if (TEXTURE_MODE() && !COLOR_MODE()) model->meshes[0].colors = cacheVertexColors;
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
    alphaMaskPunchOut(&rgba, &mask, 1);
    ImageFlipVertical(&rgba); //TODO: dont miss this, could show it properly with NDC space funkyness... still a lot of pedagogy needs to be reorganized....
    if (perspectiveCorrectTexture->id != 0)
        UpdateTexture(*perspectiveCorrectTexture, rgba.data);
    else
        *perspectiveCorrectTexture = LoadTextureFromImage(rgba);

    UnloadImage(mask);
    UnloadImage(rgba);
}

static void alphaMaskPunchOut(Image *rgba, const Image *mask, unsigned char threshold)
{
    Image maskCopy = ImageCopy(*mask);
    ImageFormat(&maskCopy, PIXELFORMAT_UNCOMPRESSED_R8G8B8A8);
    ImageFormat(rgba, PIXELFORMAT_UNCOMPRESSED_R8G8B8A8);
    unsigned char *rgbaPixels = rgba->data;
    unsigned char *maskPixels = maskCopy.data;
    int pixelCount = rgba->width * rgba->height;
    for (int i = 0; i < pixelCount; ++i)
    {
        unsigned int rchannel = maskPixels[4 * i + 0];
        unsigned int gchannel = maskPixels[4 * i + 1];
        unsigned int bchannel = maskPixels[4 * i + 2];
        unsigned int luma = (rchannel * 30u + gchannel * 59u + bchannel * 11u) / 100u;
        rgbaPixels[4 * i + 3] = (luma > threshold) ? 255 : 0;
    }
    UnloadImage(maskCopy);
}

static void fillVertexColors(Mesh *mesh)
{
    if (!mesh->colors) mesh->colors = (unsigned char *)RL_CALLOC(mesh->vertexCount, sizeof(Color));
    Color *colors = (Color *)mesh->colors;
    for (int i = 0; i < mesh->triangleCount; ++i)
    {
        const int ia = mesh->indices ? mesh->indices[3 * i + 0] : 3 * i + 0;
        const int ib = mesh->indices ? mesh->indices[3 * i + 1] : 3 * i + 1;
        const int ic = mesh->indices ? mesh->indices[3 * i + 2] : 3 * i + 2;
        colors[ia] = RED;
        colors[ib] = GREEN;
        colors[ic] = BLUE;
    }
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

    if (IsKeyDown(KEY_LEFT)) azimuth += LONG_SPEED * deltaTime;
    if (IsKeyDown(KEY_RIGHT)) azimuth -= LONG_SPEED * deltaTime;
    if (IsKeyDown(KEY_UP)) elevation += LAT_SPEED * deltaTime;
    if (IsKeyDown(KEY_DOWN)) elevation -= LAT_SPEED * deltaTime;
    if (IsKeyDown(KEY_W)) radius -= ZOOM_SPEED * deltaTime;
    if (IsKeyDown(KEY_S)) radius += ZOOM_SPEED * deltaTime;
    if (IsKeyDown(KEY_A)) rollDeltaRadians -= ROLL_SPEED * deltaTime;
    if (IsKeyDown(KEY_D)) rollDeltaRadians += ROLL_SPEED * deltaTime;
    if (IsKeyPressed(KEY_R))
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
