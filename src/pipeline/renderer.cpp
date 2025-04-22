#include "renderer.h"

#include <algorithm> //sort

#include "camera.h"
#include "../gfx/gfx.h"
#include "../gfx/shader.h"
#include "../gfx/mesh.h"
#include "../gfx/texture.h"
#include "../gfx/fbo.h"
#include "../pipeline/prefab.h"
#include "../pipeline/material.h"
#include "../pipeline/animation.h"
#include "../utils/utils.h"
#include "../extra/hdre.h"
#include "../core/ui.h"

#include "scene.h"


using namespace SCN;

//some globals
GFX::Mesh sphere;

Renderer::Renderer(const char* shader_atlas_filename)
{
	render_wireframe = false;
	render_boundaries = false;
	scene = nullptr;
	skybox_cubemap = nullptr;

	if (!GFX::Shader::LoadAtlas(shader_atlas_filename))
		exit(1);
	GFX::checkGLErrors();

	sphere.createSphere(1.0f);
	sphere.uploadToVRAM();

	Vector2ui win_wise = CORE::getWindowSize();
	shadow_FBO.setDepthOnly(1024, 1024);
	shadow_FBO1.setDepthOnly(1024, 1024);

	gbuffer.create(	win_wise.x, 
					win_wise.y,
					2,
					GL_RGBA, 
					GL_FLOAT,
					true);

	light_buffer.create(win_wise.x,
		win_wise.y,
		1,
		GL_RGBA,
		GL_FLOAT,
		true);

	quad = GFX::Mesh::getQuad();

	cone = GFX::Mesh::Get("data/meshes/cone.obj");

	sphere.createSphere(0.50f);
}

void Renderer::setupScene()
{
	if (scene->skybox_filename.size())
		skybox_cubemap = GFX::Texture::Get(std::string(scene->base_folder + "/" + scene->skybox_filename).c_str());
	else
		skybox_cubemap = nullptr;

}

void Renderer::parseNodeTree(Node* node, Camera* cam) {
	if (!node) {
		return;
	}

	// TODO: Frustrum culling

	render_calls.push_back({
		.mesh		= node->mesh,
		.model		= node->getGlobalMatrix(),
		.material	= node->material
	});

	for (Node* child_node : node->children) {
		parseNodeTree(child_node, cam);
	}
};

void Renderer::parseSceneEntities(SCN::Scene* scene, Camera* cam) {
	// HERE =====================
	// TODO: GENERATE RENDERABLES
	// ==========================

	render_calls.clear();
	scene_lights.clear();

	for (int i = 0; i < scene->entities.size(); i++) {
		BaseEntity* entity = scene->entities[i];

		if (!entity->visible) {
			continue;
		}

		if (entity->getType() == eEntityType::PREFAB) {
			PrefabEntity* ent = (PrefabEntity*)entity;
			parseNodeTree(&(ent)->root, cam);
		}
		else if (entity->getType() == eEntityType::LIGHT) {
			scene_lights.push_back((LightEntity*)entity);
		}
		// Store Prefab Entitys
		// ...
		//		Store Children Prefab Entities

		// Store Lights
		// ...
	}
	
	// Add lights to each render call
	for (int i = 0; i < render_calls.size(); i++) {
		sRenderCall* call = &render_calls[i];

		for (int j = 0; j < scene_lights.size(); j++) {
			call->lights_in_call.push_back(scene_lights[j]);
		}
	}
}

void Renderer::renderScene(SCN::Scene* scene, Camera* camera)
{
	this->scene = scene;
	setupScene();

	parseSceneEntities(scene, camera);

	//set the clear color (the background color)
	glClearColor(scene->background_color.x, scene->background_color.y, scene->background_color.z, 1.0);

	// Clear the color and the depth buffer
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	GFX::checkGLErrors();

	renderShadowMap(shadow_vp, shadow_FBO, scene_lights[3]);
	renderShadowMap(shadow_vp1, shadow_FBO1, scene_lights[0]);

	// HERE =====================
	// TODO: RENDER RENDERABLES
	// ==========================
	//gbuffer.enableAllBuffers();
	gbuffer.bind();

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	//render skybox
	if (skybox_cubemap)
		renderSkybox(skybox_cubemap);

	for (sRenderCall& render_call : render_calls) {
		renderMeshWithMaterial(Camera::current, render_call.model, render_call.mesh, render_call.material, render_call.lights_in_call);
	}

	gbuffer.unbind();

	//gbuffer.color_textures[0]->toViewport();

	// Light pass ===============

	if (use_deffered_singlepass) {
		GFX::Shader* shader = GFX::Shader::Get("deferred_light");
		if (!shader)
			return;
		shader->enable();

		// Upload light data
		vec3 light_positions[10];
		vec3 light_colors[10];
		float light_intensities[10];
		float light_types[10];

		int i = 0;
		for (LightEntity* light : scene_lights) {
			light_positions[i] = light->root.model.getTranslation();
			light_colors[i] = light->color;
			light_intensities[i] = light->intensity;
			light_types[i] = light->light_type;
			i++;
		}

		shader->setUniform("u_shadow_vp", shadow_vp);
		shader->setUniform("u_shadowmap", shadow_FBO.depth_texture, 1);

		shader->setUniform("u_shadow_vp1", shadow_vp1);
		shader->setUniform("u_shadowmap1", shadow_FBO1.depth_texture, 2);

		shader->setUniform3Array("u_light_positions", (float*)light_positions, min(10, scene_lights.size()));
		shader->setUniform3Array("u_light_colors", (float*)light_colors, min(10, scene_lights.size()));
		shader->setUniform1Array("u_light_intensities", (float*)light_intensities, min(10, scene_lights.size()));
		shader->setUniform1Array("u_light_type", (float*)light_types, min(10, scene_lights.size()));
		shader->setUniform1("u_light_count", (int)scene_lights.size());
		shader->setUniform3("u_ambient_light", vec3(0.0));

		// Upload camera uniforms
		mat4 vp_inv = camera->viewprojection_matrix;
		vp_inv.inverse();
		shader->setUniform("u_inv_vp_mat", vp_inv);
		shader->setUniform("u_camera_position", camera->eye);

		// Bing GBuffers
		shader->setTexture("u_gbuffer_albedo", gbuffer.color_textures[0], 5);
		shader->setTexture("u_gbuffer_normal", gbuffer.color_textures[1], 6);
		shader->setTexture("u_gbuffer_depth", gbuffer.depth_texture, 7);

		quad->render(GL_TRIANGLES);

		shader->disable();
	} else {
		gbuffer.depth_texture->copyTo(light_buffer.depth_texture);

		GFX::Shader* shader = GFX::Shader::Get("light_volume");
		if (!shader)
			return;
		GFX::Shader* first_shader = GFX::Shader::Get("first_pas");
		if (!first_shader)
			return;

		light_buffer.bind();

		glClear(GL_COLOR_BUFFER_BIT);

		// Upload light data
		vec3 light_positions[10];
		vec3 light_colors[10];
		float light_intensities[10];
		float light_types[10];
		vec3 light_dir[10];

		int i = 0;
		for (LightEntity* light : scene_lights) {
			light_positions[i] = light->root.model.getTranslation();
			light_colors[i] = light->color;
			light_intensities[i] = light->intensity;
			light_types[i] = light->light_type;
			light_dir[i] = light->root.model * vec3(0.0f, 0.0f, -1.0f);
			i++;
		}

		// First pass
		first_shader->enable();

		first_shader->setUniform("u_shadow_vp", shadow_vp);
		first_shader->setUniform("u_shadowmap", shadow_FBO.depth_texture, 1);
		first_shader->setUniform3("u_ambient_light", vec3(0.0));

		// Upload camera uniforms
		mat4 vp_inv = camera->viewprojection_matrix;
		vp_inv.inverse();
		first_shader->setUniform("u_inv_vp_mat", vp_inv);
		first_shader->setUniform("u_camera_position", camera->eye);

		Vector2ui win_wise = CORE::getWindowSize();
		first_shader->setUniform("u_res_inv", vec2(1.0f / win_wise.x, 1.0f / win_wise.y));
		first_shader->setUniform("u_viewprojection", camera->viewprojection_matrix);

		// Bing GBuffers
		first_shader->setTexture("u_gbuffer_albedo", gbuffer.color_textures[0], 5);
		first_shader->setTexture("u_gbuffer_normal", gbuffer.color_textures[1], 6);
		first_shader->setTexture("u_gbuffer_depth", gbuffer.depth_texture, 7);

		first_shader->setUniform("u_ambient_light", scene->ambient_light);
		first_shader->setUniform("u_light_dir", light_dir[3]);
		first_shader->setUniform("u_light_color", light_colors[3]);
		first_shader->setUniform("u_light_intensity", light_intensities[3]);

		quad->render(GL_TRIANGLES);


		first_shader->disable();

		shader->enable();

		shader->setUniform("u_shadow_vp", shadow_vp);
		shader->setUniform("u_shadowmap", shadow_FBO.depth_texture, 1);

		shader->setUniform("u_shadow_vp1", shadow_vp1);
		shader->setUniform("u_shadowmap1", shadow_FBO1.depth_texture, 2);

		shader->setUniform3Array("u_light_positions", (float*)light_positions, min(10, scene_lights.size()));
		shader->setUniform3Array("u_light_colors", (float*)light_colors, min(10, scene_lights.size()));
		shader->setUniform1Array("u_light_intensities", (float*)light_intensities, min(10, scene_lights.size()));
		shader->setUniform1Array("u_light_type", (float*)light_types, min(10, scene_lights.size()));
		shader->setUniform1("u_light_count", (int)scene_lights.size());
		
		shader->setUniform("u_inv_vp_mat", vp_inv);
		shader->setUniform("u_camera_position", camera->eye);

		shader->setUniform("u_res_inv", vec2(1.0f / win_wise.x, 1.0f / win_wise.y));
		shader->setUniform("u_viewprojection", camera->viewprojection_matrix);

		// Bing GBuffers
		shader->setTexture("u_gbuffer_albedo", gbuffer.color_textures[0], 5);
		shader->setTexture("u_gbuffer_normal", gbuffer.color_textures[1], 6);
		shader->setTexture("u_gbuffer_depth", gbuffer.depth_texture, 7);
		
		glDepthFunc(GL_GREATER);
		glDepthMask(GL_FALSE);
		glBlendFunc(GL_ONE, GL_ONE);
		glEnable(GL_BLEND);
		glFrontFace(GL_CW);
		for (int i = 0; i < scene_lights.size() -1; i++) {
			LightEntity* light = scene_lights[i];

			shader->setUniform("u_index", i);

			mat4 model;
			vec3 pos = light->root.getGlobalMatrix().getTranslation();
			model.setTranslation(pos.x, pos.y, pos.z);
			model.scale(light->max_distance, light->max_distance, light->max_distance);

			shader->setUniform("u_model", model);


			sphere.render(GL_TRIANGLES);
		}
		glDepthMask(GL_TRUE);
		glDepthFunc(GL_LEQUAL);
		glDisable(GL_BLEND);
		shader->disable();

		light_buffer.unbind();
	}

	light_buffer.color_textures[0]->toViewport();
	
	/*GFX::Shader* depth = GFX::Shader::Get("depth");
	depth->enable();
	depth->setUniform("u_camera_nearfar", vec2(scene_lights[3]->near_distance, scene_lights[3]->max_distance));
	shadow_FBO.depth_texture->toViewport(depth);
	depth->disable();*/
}

void Renderer::renderShadowMap(mat4 &shadow_viewproj, GFX::FBO& shadow_target, LightEntity* light) {

	glColorMask(false, false, false, false);

	shadow_target.bind();

	Camera light_cam;

	glClear(GL_DEPTH_BUFFER_BIT);

	float half_size = light->area / 2.0f;

	mat4 light_mat = light->root.getGlobalMatrix();
	light_cam.lookAt(light_mat.getTranslation(), light_mat * vec3(0.0f, 0.0f, -1.0f), vec3(0.0f, 1.0f, 0.0f));
	light_cam.setOrthographic(-half_size, half_size, -half_size, half_size, light->near_distance, light->max_distance);

	shadow_viewproj = light_cam.viewprojection_matrix;

	glEnable(GL_CULL_FACE);
	for (sRenderCall& render_call : render_calls) {
		renderMeshWithMaterial(&light_cam, render_call.model, render_call.mesh, render_call.material, render_call.lights_in_call);
	}
	glDisable(GL_CULL_FACE);
	glColorMask(true, true, true, true);

	shadow_target.unbind();
}


void Renderer::renderSkybox(GFX::Texture* cubemap)
{
	Camera* camera = Camera::current;

	// Apply skybox necesarry config:
	// No blending, no dpeth test, we are always rendering the skybox
	// Set the culling aproppiately, since we just want the back faces
	glDisable(GL_BLEND);
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_CULL_FACE);

	if (render_wireframe)
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

	GFX::Shader* shader = GFX::Shader::Get("skybox");
	if (!shader)
		return;
	shader->enable();

	// Center the skybox at the camera, with a big sphere
	Matrix44 m;
	m.setTranslation(camera->eye.x, camera->eye.y, camera->eye.z);
	m.scale(10, 10, 10);
	shader->setUniform("u_model", m);

	// Upload camera uniforms
	shader->setUniform("u_viewprojection", camera->viewprojection_matrix);
	shader->setUniform("u_camera_position", camera->eye);

	shader->setUniform("u_texture", cubemap, 0);

	sphere.render(GL_TRIANGLES);

	shader->disable();

	// Return opengl state to default
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glEnable(GL_DEPTH_TEST);
}

// Renders a mesh given its transform and material
void Renderer::renderMeshWithMaterial(Camera* camera, const Matrix44 model, GFX::Mesh* mesh, SCN::Material* material, std::vector<LightEntity*>& lights_to_render)
{
	//in case there is nothing to do
	if (!mesh || !mesh->getNumVertices() || !material )
		return;
    assert(glGetError() == GL_NO_ERROR);

	//define locals to simplify coding
	GFX::Shader* shader = NULL;

	glEnable(GL_DEPTH_TEST);

	//chose a shader
	shader = GFX::Shader::Get("gbuffer_fill");

    assert(glGetError() == GL_NO_ERROR);

	//no shader? then nothing to render
	if (!shader)
		return;
	shader->enable();

	material->bind(shader);

	// Upload light data
	vec3* light_positions = new vec3[lights_to_render.size()];
	vec3* light_colors = new vec3[lights_to_render.size()];
	float* light_intensities = new float[lights_to_render.size()];
	float* light_types = new float[lights_to_render.size()];

	int i = 0;
	for (LightEntity* light : lights_to_render) {
		light_positions[i] = light->root.model.getTranslation();
		light_colors[i] = light->color;
		light_intensities[i] = light->intensity;
		light_types[i] = light->light_type;
		i++;
	}

	shader->setUniform3Array("u_light_positions", (float*)light_positions, min(10, lights_to_render.size()));
	shader->setUniform3Array("u_light_colors", (float*)light_colors, min(10, lights_to_render.size()));
	shader->setUniform1Array("u_light_intensities", (float*)light_intensities, min(10, lights_to_render.size()));
	shader->setUniform1Array("u_light_type", (float*)light_types, min(10, lights_to_render.size()));
	shader->setUniform1("u_light_count", (int) lights_to_render.size());
	shader->setUniform3("u_ambient_light", vec3(0.0));

	shader->setUniform("u_shadow_vp", shadow_vp);
	shader->setUniform("u_shadowmap", shadow_FBO.depth_texture, 1);

	shader->setUniform("u_shadow_vp1", shadow_vp1);
	shader->setUniform("u_shadowmap1", shadow_FBO1.depth_texture, 2);

	//upload uniforms
	shader->setUniform("u_model", model);

	// Upload camera uniforms
	shader->setUniform("u_viewprojection", camera->viewprojection_matrix);
	shader->setUniform("u_camera_position", camera->eye);

	// Upload time, for cool shader effects
	float t = getTime();
	shader->setUniform("u_time", t );

	// Render just the verticies as a wireframe
	if (render_wireframe)
		glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );

	//do the draw call that renders the mesh into the screen
	mesh->render(GL_TRIANGLES);

	//disable shader
	shader->disable();

	//set the render state as it was before to avoid problems with future renders
	glDisable(GL_BLEND);
	glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );

	/*delete[] light_positions;
	delete[] light_intensities;
	delete[] light_types;
	delete[] light_colors;*/
}

#ifndef SKIP_IMGUI

void Renderer::showUI()
{
		
	ImGui::Checkbox("Wireframe", &render_wireframe);
	ImGui::Checkbox("Boundaries", &render_boundaries);

	//add here your stuff
	//...
}

#else
void Renderer::showUI() {}
#endif