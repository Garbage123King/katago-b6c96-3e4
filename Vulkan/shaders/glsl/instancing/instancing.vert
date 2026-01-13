#version 450

layout (location = 0) in vec3 inPos;
layout (location = 1) in vec3 inNormal;
layout (location = 2) in vec2 inUV;
layout (location = 3) in vec3 inColor;

// 实例属性
layout (location = 4) in vec3 instancePos;
layout (location = 5) in vec3 instanceRot;
layout (location = 6) in float instanceScale;
layout (location = 7) in float instanceColorValue; 

layout (binding = 0) uniform UBO {
    mat4 projection;
    mat4 view;
} ubo;

layout (location = 0) out float outColor; 

void main() {
    outColor = instanceColorValue;

    // 直接在局部坐标上应用缩放
    vec3 localPos = inPos * (instanceScale == 0.0 ? 1.0 : instanceScale);
    
    // 如果不需要旋转，直接将缩放后的坐标加上实例世界坐标
    vec3 worldPos = localPos + instancePos;

    gl_Position = ubo.projection * ubo.view * vec4(worldPos, 1.0);
}