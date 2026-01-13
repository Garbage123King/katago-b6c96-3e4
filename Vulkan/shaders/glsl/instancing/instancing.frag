#version 450

layout (location = 0) in float inColor;
layout (location = 0) out vec4 outFragColor;

void main() {
    // 输出灰度颜色：随机数 0 为黑，1 为白
    outFragColor = vec4(vec3(inColor), 1.0);
}