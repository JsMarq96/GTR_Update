#pragma once
#include "scene.h"
#include "prefab.h"
#include "../gfx/fbo.h"
#include "light.h"

//forward declarations
class Camera;
class Skeleton;
namespace GFX {
	class Shader;
	class Mesh;
	class FBO;
}

namespace SCN {

	class Prefab;
	class Material;

	// This class is in charge of rendering anything in our system.
	// Separating the render from anything else makes the code cleaner
	class Renderer
	{
	public:

		struct sRenderCall {
			GFX::Mesh		*mesh = nullptr;
			Matrix44		model;
			Material		*material = nullptr;

			std::vector<LightEntity*> lights_in_call;
		};

		std::vector<LightEntity*> scene_lights;
		std::vector<sRenderCall> render_calls;

		GFX::Mesh* quad;

		GFX::FBO shadow_FBO;
		mat4 shadow_vp;

		GFX::FBO shadow_FBO1;
		mat4 shadow_vp1;

		bool render_wireframe;
		bool render_boundaries;

		GFX::FBO gbuffer;
		GFX::FBO light_buffer;

		bool use_deffered_singlepass = false;

		GFX::Texture* skybox_cubemap;

		SCN::Scene* scene;
		GFX::Mesh* cone;
		GFX::Mesh sphere;

		//updated every frame
		Renderer(const char* shaders_atlas_filename );

		//just to be sure we have everything ready for the rendering
		void setupScene();

		//add here your functions
		//...
		void parseNodeTree(Node* node, Camera* cam);

		void parseSceneEntities(SCN::Scene* scene, Camera* camera);

		//renders several elements of the scene
		void renderScene(SCN::Scene* scene, Camera* camera);

		//render the skybox
		void renderSkybox(GFX::Texture* cubemap);

		void renderShadowMap(mat4& shadow_viewproj, GFX::FBO& shadow_target, LightEntity* light);

		//to render one mesh given its material and transformation matrix
		void renderMeshWithMaterial(Camera* camera, const Matrix44 model, GFX::Mesh* mesh, SCN::Material* material, std::vector<LightEntity*> &lights_to_render);

		void showUI();
	};

};