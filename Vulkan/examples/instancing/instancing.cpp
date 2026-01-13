#include "vulkanexamplebase.h"
#include <random>

// 数组维度定义
#define DIM_W 100
#define DIM_X 96
#define DIM_Y 19
#define DIM_Z 19
#define INSTANCE_COUNT (DIM_W * DIM_X * DIM_Y * DIM_Z)

class VulkanExample : public VulkanExampleBase
{
public:
	struct CubeVertex {
		glm::vec3 pos;
		glm::vec3 normal;
		glm::vec2 uv;
		glm::vec3 color;
	};

	// 修正后的结构体定义
	struct MeshBuffer {
		VkBuffer buffer = VK_NULL_HANDLE;
		VkDeviceMemory memory = VK_NULL_HANDLE;
		uint32_t count = 0;
	};

	struct CubeModel {
		MeshBuffer vertices;
		MeshBuffer indices;
	} cubeModel;

	struct InstanceData {
		glm::vec3 pos;
		glm::vec3 rot;
		float scale;
		float colorValue;
	};

	struct {
		VkBuffer buffer = VK_NULL_HANDLE;
		VkDeviceMemory memory = VK_NULL_HANDLE;
		size_t size = 0;
	} instanceBuffer;

	struct UniformData {
		glm::mat4 projection;
		glm::mat4 view;
	} uniformData;
	std::array<vks::Buffer, maxConcurrentFrames> uniformBuffers;

	VkPipelineLayout pipelineLayout = VK_NULL_HANDLE;
	struct {
		VkPipeline instancedCubes = VK_NULL_HANDLE;
	} pipelines;

	VkDescriptorSetLayout descriptorSetLayout = VK_NULL_HANDLE;
	std::array<VkDescriptorSet, maxConcurrentFrames> descriptorSets;

	VulkanExample() : VulkanExampleBase()
	{
		title = "Vulkan Instance Array [100][96][19][19]";
		camera.type = Camera::CameraType::firstperson;
		//camera.setPosition(glm::vec3(0.0f, 0.0f, -50.0f));
		//camera.setRotation(glm::vec3(0.0f));
		//camera.setPerspective(60.0f, (float)width / (float)height, 0.1f, 1000.0f);
		// 把相机放到方块阵列的侧前方
		// X 设为阵列中心 (96*3/2 = 144)
		// Y 设为阵列中心 (19*3/2 = 28)
		// Z 设为负数，这样你就能从外面往里看
		// 默认速度通常是 1.0f 左右，根据你的场景大小，建议提升 50-100 倍
		camera.movementSpeed = 50.0f;
		camera.setPosition(glm::vec3(144.0f, 28.0f, -150.0f));
		camera.setRotation(glm::vec3(0.0f, 0.0f, 0.0f));
		camera.setPerspective(60.0f, (float)width / (float)height, 0.1f, 5000.0f); // 远裁剪面加大到 5000
	}

	~VulkanExample()
	{
		if (device) {
			vkDestroyPipeline(device, pipelines.instancedCubes, nullptr);
			vkDestroyPipelineLayout(device, pipelineLayout, nullptr);
			vkDestroyDescriptorSetLayout(device, descriptorSetLayout, nullptr);
			vkDestroyBuffer(device, instanceBuffer.buffer, nullptr);
			vkFreeMemory(device, instanceBuffer.memory, nullptr);
			vkDestroyBuffer(device, cubeModel.vertices.buffer, nullptr);
			vkFreeMemory(device, cubeModel.vertices.memory, nullptr);
			vkDestroyBuffer(device, cubeModel.indices.buffer, nullptr);
			vkFreeMemory(device, cubeModel.indices.memory, nullptr);
			for (auto& buffer : uniformBuffers) {
				buffer.destroy();
			}
		}
	}

	void createCubeModel() {
		// 标准立方体数据
		std::vector<CubeVertex> vertices = {
			{ {-0.5f, -0.5f,  0.5f}, {0,0,1}, {0,0}, {1,1,1} }, { { 0.5f, -0.5f,  0.5f}, {0,0,1}, {1,0}, {1,1,1} },
			{ { 0.5f,  0.5f,  0.5f}, {0,0,1}, {1,1}, {1,1,1} }, { {-0.5f,  0.5f,  0.5f}, {0,0,1}, {0,1}, {1,1,1} },
			{ {-0.5f, -0.5f, -0.5f}, {0,0,-1}, {1,0}, {1,1,1} }, { { 0.5f, -0.5f, -0.5f}, {0,0,-1}, {0,0}, {1,1,1} },
			{ { 0.5f,  0.5f, -0.5f}, {0,0,-1}, {0,1}, {1,1,1} }, { {-0.5f,  0.5f, -0.5f}, {0,0,-1}, {1,1}, {1,1,1} },
			// ... (为了简洁，这里逻辑上应包含所有6个面， Sascha框架通常有辅助函数，这里手动简化)
		};
		// 补全所有面数据以保证渲染正确
		vertices = {
			{{-1,-1, 1},{0,0, 1},{0,0},{1,1,1}}, {{ 1,-1, 1},{0,0, 1},{1,0},{1,1,1}}, {{ 1, 1, 1},{0,0, 1},{1,1},{1,1,1}}, {{-1, 1, 1},{0,0, 1},{0,1},{1,1,1}},
			{{-1,-1,-1},{0,0,-1},{1,0},{1,1,1}}, {{ 1,-1,-1},{0,0,-1},{0,0},{1,1,1}}, {{ 1, 1,-1},{0,0,-1},{0,1},{1,1,1}}, {{-1, 1,-1},{0,0,-1},{1,1},{1,1,1}},
			{{ 1,-1, 1},{1,0, 0},{0,0},{1,1,1}}, {{ 1,-1,-1},{1,0, 0},{1,0},{1,1,1}}, {{ 1, 1,-1},{1,0, 0},{1,1},{1,1,1}}, {{ 1, 1, 1},{1,0, 0},{0,1},{1,1,1}},
			{{-1,-1, 1},{-1,0,0},{1,0},{1,1,1}}, {{-1,-1,-1},{-1,0,0},{0,0},{1,1,1}}, {{-1, 1,-1},{-1,0,0},{0,1},{1,1,1}}, {{-1, 1, 1},{-1,0,0},{1,1},{1,1,1}},
			{{-1, 1, 1},{0, 1,0},{0,0},{1,1,1}}, {{ 1, 1, 1},{0, 1,0},{1,0},{1,1,1}}, {{ 1, 1,-1},{0, 1,0},{1,1},{1,1,1}}, {{-1, 1,-1},{0, 1,0},{0,1},{1,1,1}},
			{{-1,-1, 1},{0,-1,0},{0,1},{1,1,1}}, {{ 1,-1, 1},{0,-1,0},{1,1},{1,1,1}}, {{ 1,-1,-1},{0,-1,0},{1,0},{1,1,1}}, {{-1,-1,-1},{0,-1,0},{0,0},{1,1,1}}
		};

		std::vector<uint32_t> indices = {
			0,1,2, 2,3,0, 4,5,6, 6,7,4, 8,9,10, 10,11,8,
			12,13,14, 14,15,12, 16,17,18, 18,19,16, 20,21,22, 22,23,20
		};

		cubeModel.vertices.count = (uint32_t)vertices.size();
		cubeModel.indices.count = (uint32_t)indices.size();

		// 使用框架封装的创建缓冲函数 (注意参数匹配)
		VK_CHECK_RESULT(vulkanDevice->createBuffer(
			VK_BUFFER_USAGE_VERTEX_BUFFER_BIT,
			VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
			cubeModel.vertices.count * sizeof(CubeVertex),
			&cubeModel.vertices.buffer,
			&cubeModel.vertices.memory,
			vertices.data()));

		VK_CHECK_RESULT(vulkanDevice->createBuffer(
			VK_BUFFER_USAGE_INDEX_BUFFER_BIT,
			VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
			cubeModel.indices.count * sizeof(uint32_t),
			&cubeModel.indices.buffer,
			&cubeModel.indices.memory,
			indices.data()));
	}

	void prepareInstanceData()
	{
		std::vector<InstanceData> instanceData(INSTANCE_COUNT);
		std::default_random_engine rndGenerator((unsigned)time(NULL));
		std::uniform_real_distribution<float> uniformDist(0.0, 1.0);

		float spacing = 3.0f;
		for (int w = 0; w < DIM_W; w++) {
			for (int x = 0; x < DIM_X; x++) {
				for (int y = 0; y < DIM_Y; y++) {
					for (int z = 0; z < DIM_Z; z++) {
						uint32_t idx = w * (DIM_X * DIM_Y * DIM_Z) + x * (DIM_Y * DIM_Z) + y * DIM_Z + z;
						instanceData[idx].pos = glm::vec3(x * spacing, y * spacing, z * spacing + (w * 100.0f));
						instanceData[idx].rot = glm::vec3(0.0f);
						instanceData[idx].scale = 1.0f;
						instanceData[idx].colorValue = uniformDist(rndGenerator);
					}
				}
			}
		}

		instanceBuffer.size = instanceData.size() * sizeof(InstanceData);
		// 创建 Device Local 缓冲区（略去 Staging 步骤以简化，直接创建 Host Visible）
		VK_CHECK_RESULT(vulkanDevice->createBuffer(
			VK_BUFFER_USAGE_VERTEX_BUFFER_BIT,
			VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT,
			instanceBuffer.size,
			&instanceBuffer.buffer,
			&instanceBuffer.memory,
			instanceData.data()));
	}

	void setupDescriptors()
	{
		std::vector<VkDescriptorPoolSize> poolSizes = { vks::initializers::descriptorPoolSize(VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, maxConcurrentFrames) };
		VkDescriptorPoolCreateInfo descriptorPoolInfo = vks::initializers::descriptorPoolCreateInfo(poolSizes, maxConcurrentFrames);
		VK_CHECK_RESULT(vkCreateDescriptorPool(device, &descriptorPoolInfo, nullptr, &descriptorPool));

		std::vector<VkDescriptorSetLayoutBinding> setLayoutBindings = { vks::initializers::descriptorSetLayoutBinding(VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, VK_SHADER_STAGE_VERTEX_BIT, 0) };
		VkDescriptorSetLayoutCreateInfo descriptorLayout = vks::initializers::descriptorSetLayoutCreateInfo(setLayoutBindings);
		VK_CHECK_RESULT(vkCreateDescriptorSetLayout(device, &descriptorLayout, nullptr, &descriptorSetLayout));

		for (uint32_t i = 0; i < maxConcurrentFrames; i++) {
			VkDescriptorSetAllocateInfo allocInfo = vks::initializers::descriptorSetAllocateInfo(descriptorPool, &descriptorSetLayout, 1);
			VK_CHECK_RESULT(vkAllocateDescriptorSets(device, &allocInfo, &descriptorSets[i]));
			VkWriteDescriptorSet writeDescriptorSet = vks::initializers::writeDescriptorSet(descriptorSets[i], VK_DESCRIPTOR_TYPE_UNIFORM_BUFFER, 0, &uniformBuffers[i].descriptor);
			vkUpdateDescriptorSets(device, 1, &writeDescriptorSet, 0, nullptr);
		}
	}

	void preparePipelines()
	{
		VkPipelineLayoutCreateInfo pipelineLayoutCI = vks::initializers::pipelineLayoutCreateInfo(&descriptorSetLayout, 1);
		VK_CHECK_RESULT(vkCreatePipelineLayout(device, &pipelineLayoutCI, nullptr, &pipelineLayout));

		VkPipelineInputAssemblyStateCreateInfo inputAssemblyState = vks::initializers::pipelineInputAssemblyStateCreateInfo(VK_PRIMITIVE_TOPOLOGY_TRIANGLE_LIST, 0, VK_FALSE);
		VkPipelineRasterizationStateCreateInfo rasterizationState = vks::initializers::pipelineRasterizationStateCreateInfo(VK_POLYGON_MODE_FILL, VK_CULL_MODE_BACK_BIT, VK_FRONT_FACE_COUNTER_CLOCKWISE, 0);
		VkPipelineColorBlendAttachmentState blendAttachmentState = vks::initializers::pipelineColorBlendAttachmentState(0xf, VK_FALSE);
		VkPipelineColorBlendStateCreateInfo colorBlendState = vks::initializers::pipelineColorBlendStateCreateInfo(1, &blendAttachmentState);
		VkPipelineDepthStencilStateCreateInfo depthStencilState = vks::initializers::pipelineDepthStencilStateCreateInfo(VK_TRUE, VK_TRUE, VK_COMPARE_OP_LESS_OR_EQUAL);
		VkPipelineViewportStateCreateInfo viewportState = vks::initializers::pipelineViewportStateCreateInfo(1, 1, 0);
		VkPipelineMultisampleStateCreateInfo multisampleState = vks::initializers::pipelineMultisampleStateCreateInfo(VK_SAMPLE_COUNT_1_BIT, 0);
		std::vector<VkDynamicState> dynamicStateEnables = { VK_DYNAMIC_STATE_VIEWPORT, VK_DYNAMIC_STATE_SCISSOR };
		VkPipelineDynamicStateCreateInfo dynamicState = vks::initializers::pipelineDynamicStateCreateInfo(dynamicStateEnables);

		// 顶点输入属性
		std::vector<VkVertexInputBindingDescription> bindingDescriptions = {
			vks::initializers::vertexInputBindingDescription(0, sizeof(CubeVertex), VK_VERTEX_INPUT_RATE_VERTEX),
			vks::initializers::vertexInputBindingDescription(1, sizeof(InstanceData), VK_VERTEX_INPUT_RATE_INSTANCE)
		};
		std::vector<VkVertexInputAttributeDescription> attributeDescriptions = {
			vks::initializers::vertexInputAttributeDescription(0, 0, VK_FORMAT_R32G32B32_SFLOAT, offsetof(CubeVertex, pos)),
			vks::initializers::vertexInputAttributeDescription(1, 4, VK_FORMAT_R32G32B32_SFLOAT, offsetof(InstanceData, pos)),
			vks::initializers::vertexInputAttributeDescription(1, 7, VK_FORMAT_R32_SFLOAT, offsetof(InstanceData, colorValue))
		};
		VkPipelineVertexInputStateCreateInfo inputState = vks::initializers::pipelineVertexInputStateCreateInfo();
		inputState.vertexBindingDescriptionCount = (uint32_t)bindingDescriptions.size();
		inputState.pVertexBindingDescriptions = bindingDescriptions.data();
		inputState.vertexAttributeDescriptionCount = (uint32_t)attributeDescriptions.size();
		inputState.pVertexAttributeDescriptions = attributeDescriptions.data();

		VkGraphicsPipelineCreateInfo pipelineCI = vks::initializers::pipelineCreateInfo(pipelineLayout, renderPass);
		pipelineCI.pVertexInputState = &inputState;
		pipelineCI.pInputAssemblyState = &inputAssemblyState;
		pipelineCI.pRasterizationState = &rasterizationState;
		pipelineCI.pColorBlendState = &colorBlendState;
		pipelineCI.pMultisampleState = &multisampleState;
		pipelineCI.pViewportState = &viewportState;
		pipelineCI.pDepthStencilState = &depthStencilState;
		pipelineCI.pDynamicState = &dynamicState;

		std::array<VkPipelineShaderStageCreateInfo, 2> shaderStages;
		shaderStages[0] = loadShader(getShadersPath() + "instancing/instancing.vert.spv", VK_SHADER_STAGE_VERTEX_BIT);
		shaderStages[1] = loadShader(getShadersPath() + "instancing/instancing.frag.spv", VK_SHADER_STAGE_FRAGMENT_BIT);
		pipelineCI.stageCount = 2;
		pipelineCI.pStages = shaderStages.data();

		VK_CHECK_RESULT(vkCreateGraphicsPipelines(device, pipelineCache, 1, &pipelineCI, nullptr, &pipelines.instancedCubes));
	}

	void buildCommandBuffer()
	{
		VkCommandBufferBeginInfo cmdBufInfo = vks::initializers::commandBufferBeginInfo();
		VkClearValue clearValues[2];
		clearValues[0].color = { { 0.1f, 0.1f, 0.1f, 1.0f } };
		clearValues[1].depthStencil = { 1.0f, 0 };

		VkRenderPassBeginInfo renderPassBeginInfo = vks::initializers::renderPassBeginInfo();
		renderPassBeginInfo.renderPass = renderPass;
		renderPassBeginInfo.renderArea.extent.width = width;
		renderPassBeginInfo.renderArea.extent.height = height;
		renderPassBeginInfo.clearValueCount = 2;
		renderPassBeginInfo.pClearValues = clearValues;
		renderPassBeginInfo.framebuffer = frameBuffers[currentImageIndex];

		VK_CHECK_RESULT(vkBeginCommandBuffer(drawCmdBuffers[currentBuffer], &cmdBufInfo));
		vkCmdBeginRenderPass(drawCmdBuffers[currentBuffer], &renderPassBeginInfo, VK_SUBPASS_CONTENTS_INLINE);

		VkViewport viewport = vks::initializers::viewport((float)width, (float)height, 0.0f, 1.0f);
		vkCmdSetViewport(drawCmdBuffers[currentBuffer], 0, 1, &viewport);
		VkRect2D scissor = vks::initializers::rect2D(width, height, 0, 0);
		vkCmdSetScissor(drawCmdBuffers[currentBuffer], 0, 1, &scissor);

		vkCmdBindDescriptorSets(drawCmdBuffers[currentBuffer], VK_PIPELINE_BIND_POINT_GRAPHICS, pipelineLayout, 0, 1, &descriptorSets[currentBuffer], 0, nullptr);
		vkCmdBindPipeline(drawCmdBuffers[currentBuffer], VK_PIPELINE_BIND_POINT_GRAPHICS, pipelines.instancedCubes);

		VkDeviceSize offsets[1] = { 0 };
		// 修正这里的绑定方式
		vkCmdBindVertexBuffers(drawCmdBuffers[currentBuffer], 0, 1, &cubeModel.vertices.buffer, offsets);
		vkCmdBindVertexBuffers(drawCmdBuffers[currentBuffer], 1, 1, &instanceBuffer.buffer, offsets);
		vkCmdBindIndexBuffer(drawCmdBuffers[currentBuffer], cubeModel.indices.buffer, 0, VK_INDEX_TYPE_UINT32);

		// vkCmdDrawIndexed 确实是 5 个参数：(cmd, indexCount, instanceCount, firstIndex, vertexOffset, firstInstance)
		vkCmdDrawIndexed(drawCmdBuffers[currentBuffer], cubeModel.indices.count, INSTANCE_COUNT, 0, 0, 0);

		drawUI(drawCmdBuffers[currentBuffer]);
		vkCmdEndRenderPass(drawCmdBuffers[currentBuffer]);
		VK_CHECK_RESULT(vkEndCommandBuffer(drawCmdBuffers[currentBuffer]));
	}

	void prepareUniformBuffers()
	{
		for (auto& buffer : uniformBuffers) {
			VK_CHECK_RESULT(vulkanDevice->createBuffer(VK_BUFFER_USAGE_UNIFORM_BUFFER_BIT, VK_MEMORY_PROPERTY_HOST_VISIBLE_BIT | VK_MEMORY_PROPERTY_HOST_COHERENT_BIT, &buffer, sizeof(UniformData)));
			VK_CHECK_RESULT(buffer.map());
		}
	}

	void updateUniformBuffers()
	{
		uniformData.projection = camera.matrices.perspective;
		uniformData.view = camera.matrices.view;
		memcpy(uniformBuffers[currentBuffer].mapped, &uniformData, sizeof(uniformData));
	}

	void prepare()
	{
		VulkanExampleBase::prepare();
		createCubeModel();
		prepareInstanceData();
		prepareUniformBuffers();
		setupDescriptors();
		preparePipelines();
		prepared = true;
	}

	virtual void render()
	{
		if (!prepared) return;
		VulkanExampleBase::prepareFrame();
		updateUniformBuffers();
		buildCommandBuffer();
		VulkanExampleBase::submitFrame();
		//std::cout << "Camera Pos: " << camera.position.z << std::endl;
	}
};

VULKAN_EXAMPLE_MAIN()