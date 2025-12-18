// raylib_example.c - 完整实现你的所有需求
#include "raylib.h"

int main() {
    // 初始化窗口 - 一行代码
    InitWindow(800, 600, "Simple3dQube");
    
    // 设置相机 - 几行代码
    Camera camera = { 0 };
    camera.position = (Vector3){ 10.0f, 10.0f, 10.0f };
    camera.target = (Vector3){ 0.0f, 0.0f, 0.0f };
    camera.up = (Vector3){ 0.0f, 1.0f, 0.0f };
    camera.fovy = 45.0f;
    camera.projection = CAMERA_PERSPECTIVE;
    
    // 创建3D方块模型
    Vector3 cubePosition = { 0.0f, 0.0f, 0.0f };
    
    // 加载字体（内置默认字体）
    Font font = GetFontDefault();
    
    // 动画变量
    float rotation = 0.0f;
    
    // 设置帧率
    SetTargetFPS(60);
    
    // 主循环
    while (!WindowShouldClose()) {
        // 更新动画
        rotation += 0.5f;
        
        // 相机控制（按WASD移动）
        UpdateCamera(&camera, CAMERA_ORBITAL);
        
        // 开始绘制
        BeginDrawing();
            ClearBackground(RAYWHITE);
            
            // 3D绘制
            BeginMode3D(camera);
                // 绘制网格（辅助参考）
                DrawGrid(10, 1.0f);
                
                // 绘制3D方块 - 带光照
                DrawCube(cubePosition, 2.0f, 2.0f, 2.0f, RED);
                DrawCubeWires(cubePosition, 2.0f, 2.0f, 2.0f, MAROON);
                
                // 绘制旋转的方块
                DrawCube((Vector3){4.0f, 0.0f, 0.0f}, 1.0f, 1.0f, 1.0f, BLUE);
                DrawCubeWires((Vector3){4.0f, 0.0f, 0.0f}, 1.0f, 1.0f, 1.0f, DARKBLUE);
                
                // 绘制金字塔
                DrawCube((Vector3){-4.0f, 0.0f, 0.0f}, 1.5f, 3.0f, 1.5f, GREEN);
                DrawCubeWires((Vector3){-4.0f, 0.0f, 0.0f}, 1.5f, 3.0f, 1.5f, DARKGREEN);
                
            EndMode3D();
            
            // 2D文字显示
            DrawText("Simple 3D Demo", 10, 10, 20, DARKGRAY);
            DrawText("WASD: move camera  mouse: rotate vision", 10, 40, 20, DARKGRAY);
            DrawText(TextFormat("Qube rotate angle: %.1f", rotation), 10, 70, 20, DARKGRAY);
            DrawFPS(10, 100);
            
        EndDrawing();
    }
    
    // 清理资源
    CloseWindow();
    
    return 0;
}