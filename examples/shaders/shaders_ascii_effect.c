/*******************************************************************************************
*
*   raylib [shaders] example - ascii effect
*
*   Example complexity rating: [★☆☆☆] 1/4
*
*   Example originally created with raylib 5.5, last time updated with raylib 5.6
*
*   Example contributed by Maicon Santana (@maiconpintoabreu) and reviewed by Ramon Santamaria (@raysan5)
*
*   Example licensed under an unmodified zlib/libpng license, which is an OSI-certified,
*   BSD-like license that allows static linking with closed source software
*
*   Copyright (c) 2025-2025 Maicon Santana (@maiconpintoabreu)
*
********************************************************************************************/

#include "raylib.h"

#if defined(PLATFORM_DESKTOP)
    #define GLSL_VERSION            330
#else   // PLATFORM_ANDROID, PLATFORM_WEB
    #define GLSL_VERSION            100
#endif

//------------------------------------------------------------------------------------
// Program main entry point
//------------------------------------------------------------------------------------
int main(void)
{
    // Initialization
    //--------------------------------------------------------------------------------------
    const int screenWidth = 800;
    const int screenHeight = 450;

    InitWindow(screenWidth, screenHeight, "raylib [shaders] example - ascii effect");

    // Texture to test static drawing
    Texture2D fudesumi = LoadTexture("resources/fudesumi.png");
    // Texture to test moving drawing
    Texture2D raysan = LoadTexture("resources/raysan.png");

    // Load shader to be used on postprocessing
    Shader shader = LoadShader(0, TextFormat("resources/shaders/glsl%i/ascii.fs", GLSL_VERSION));

    // These locations are used to send data to the GPU.
    int resolutionLoc = GetShaderLocation(shader, "resolution");
    int fontSizeLoc = GetShaderLocation(shader, "fontSize");

    // Set the character size for the ASCII effect
    float fontSize = 4.0f;

    // Send the updated values to the shader
    float resolution[2] = { (float)screenWidth, (float)screenHeight };
    SetShaderValue(shader, resolutionLoc, resolution, SHADER_UNIFORM_VEC2);
    SetShaderValue(shader, fontSizeLoc, &fontSize, SHADER_UNIFORM_FLOAT);

    Vector2 circlePos = (Vector2){40.0f, (float)screenHeight * 0.5f};
    float circleSpeed = 1.0f;

    // RenderTexture to apply the postprocessing later
    RenderTexture2D target = LoadRenderTexture(screenWidth, screenHeight);

    SetTargetFPS(60);                   // Set our game to run at 60 frames-per-second
    //--------------------------------------------------------------------------------------

    // Main game loop
    while (!WindowShouldClose())    // Detect window close button or ESC key
    {
        // Update
        //----------------------------------------------------------------------------------
        if (circlePos.x > 200.0f || circlePos.x < 40.0f) {
            circleSpeed *= -1;
        }
        circlePos.x += circleSpeed;

        // Draw
        //----------------------------------------------------------------------------------

        BeginTextureMode(target);
            ClearBackground(WHITE); // The background of the scene itself
            
            DrawTexture(fudesumi, 500, -30, WHITE);    // Using custom shader
            DrawTextureV(raysan, circlePos, WHITE);

        EndTextureMode();
        BeginDrawing();

            ClearBackground(RAYWHITE);

            BeginShaderMode(shader);

                // Draw the scene texture (that we rendered earlier) to the screen.
                // The shader will process every pixel of this texture.
                DrawTextureRec(target.texture, 
                               (Rectangle){ 0, 0, (float)target.texture.width, (float)-target.texture.height }, 
                               (Vector2){ 0, 0 }, 
                               WHITE);
            EndShaderMode();

            DrawRectangle(0, 0, screenWidth, 40, BLACK);
            DrawText("Ascii effect", 120, 10, 20, LIGHTGRAY);
            DrawFPS(10, 10);

        EndDrawing();
        //----------------------------------------------------------------------------------
    }

    // De-Initialization
    //--------------------------------------------------------------------------------------
    UnloadShader(shader);       // Unload shader
    UnloadTexture(fudesumi);    // Unload texture
    UnloadTexture(raysan);      // Unload texture

    CloseWindow();        // Close window and OpenGL context
    //--------------------------------------------------------------------------------------

    return 0;
}
