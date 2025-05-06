//example of some shaders compiled
flat basic.vs flat.fs
texture basic.vs texture.fs
plain basic.vs plain.fs
skybox basic.vs skybox.fs
depth quad.vs depth.fs
deferred_light quad.vs deferred_light.fs
multi basic.vs multi.fs
gbuffer_fill basic.vs gbuffer_fill.fs
light_volume basic.vs light_volume.fs
first_pas quad.vs first_pass.fs
ssao_pass quad.vs ssao_pass.fs


phong basic.vs phong.fs

\PBR

float PI = 3.14159265;

vec3 F_Schlick(vec3 v, vec3 n, vec3 f0) {
	return f0 + (1.0 - f0) * pow(1.0 -  clamp(dot(n, v), 0.0001, 1.0), 5.0);
}

float D_ggx(float alpha, vec3 n, vec3 h) {
	float alpha_2 = alpha * alpha;
	float n_dot_h = clamp(dot(n,h), 0.000, 1.0);
	float dem = (n_dot_h*n_dot_h) * (alpha_2 - 1.0) + 1.0;
	return alpha_2 / (PI * dem * dem);
}

float schlick_ggx(float k, vec3 v, vec3 n) {
	float n_dot_v = clamp(dot(n,v), 0.0001, 1.0);
	return n_dot_v / (n_dot_v * (1.0-k) + k);
}

float G_smith(float alpha, vec3 l, vec3 v, vec3 n) {
	float k = alpha / 2.0;
	return schlick_ggx(k, l, n) * schlick_ggx(k, v, n);
}

vec3 brdf_cook_torrance(vec3 albedo, float roughness, float metalness, vec3 n, vec3 h, vec3 l, vec3 v) {
	vec3 F0 = mix(vec3(0.04), albedo, metalness);
	float alpha = roughness * roughness; 

	vec3 diffuse = albedo / PI;
	
	vec3 F = F_Schlick(h, v, F0);
	float G = G_smith(alpha, l, v, n);
	float D = D_ggx(alpha, n, h);

	vec3 specular = (F * G * D) / ((4.0 * clamp(dot(n, l), 0.0001, 1.0) * clamp(dot(n, v), 0.0001, 1.0)));

	return diffuse + specular; 
}

\basic.vs

#version 330 core

in vec3 a_vertex;
in vec3 a_normal;
in vec2 a_coord;
in vec4 a_color;

uniform vec3 u_camera_pos;

uniform mat4 u_model;
uniform mat4 u_viewprojection;

//this will store the color for the pixel shader
out vec3 v_position;
out vec3 v_world_position;
out vec3 v_normal;
out vec2 v_uv;
out vec4 v_color;

uniform float u_time;

void main()
{	
	//calcule the normal in camera space (the NormalMatrix is like ViewMatrix but without traslation)
	v_normal = (u_model * vec4( a_normal, 0.0) ).xyz;
	
	//calcule the vertex in object space
	v_position = a_vertex;
	v_world_position = (u_model * vec4( v_position, 1.0) ).xyz;
	
	//store the color in the varying var to use it from the pixel shader
	v_color = a_color;

	//store the texture coordinates
	v_uv = a_coord;

	//calcule the position of the vertex using the matrices
	gl_Position = u_viewprojection * vec4( v_world_position, 1.0 );
}

\quad.vs

#version 330 core

in vec3 a_vertex;
in vec2 a_coord;
out vec2 v_uv;

void main()
{	
	v_uv = a_coord;
	gl_Position = vec4( a_vertex, 1.0 );
}


\flat.fs

#version 330 core

uniform vec4 u_color;

out vec4 FragColor;

void main()
{
	FragColor = u_color;
}


\plain.fs

#version 330 core

out vec4 FragColor;

void main()
{
	FragColor = vec4(0.0, 0.0, 0.0, 1.0);
}

\texture.fs

#version 330 core

in vec3 v_position;
in vec3 v_world_position;
in vec3 v_normal;
in vec2 v_uv;
in vec4 v_color;

uniform vec4 u_color;
uniform sampler2D u_texture;
uniform float u_time;
uniform float u_alpha_cutoff;

out vec4 FragColor;

void main()
{
	vec2 uv = v_uv;
	vec4 color = u_color;
	color *= texture2D( u_texture, uv );

	if(color.a < u_alpha_cutoff)
		discard;

	FragColor = color;
}

\gbuffer_fill.fs

#version 330 core

in vec3 v_position;
in vec3 v_world_position;
in vec3 v_normal;
in vec2 v_uv;
in vec4 v_color;

uniform vec4 u_color;
uniform sampler2D u_texture;
uniform sampler2D u_normalmap;
uniform sampler2D u_mrmap;
uniform float u_time;
uniform float u_alpha_cutoff;
uniform vec3 u_camera_position;

layout(location = 0) out vec4 gbuffer_albedo;
layout(location = 1) out vec4 gbuffer_normal_mat;

void main()
{
	vec2 uv = v_uv;
	// Base color
	vec4 color = u_color;
	color *= texture2D( u_texture, uv );

	vec4 material = texture2D(u_mrmap, uv);

	if(color.a < u_alpha_cutoff)
		discard;

	vec3 N = normalize(v_normal);

	gbuffer_albedo = vec4(color.xyz, material.g); // Roughness
	gbuffer_normal_mat = vec4(N * 0.5 + 0.5, material.b); // Metalness
}

\deferred_light.fs
#version 330 core

in vec2 v_uv;

uniform mat4 u_inv_vp_mat;

uniform sampler2D u_gbuffer_albedo;
uniform sampler2D u_gbuffer_normal;
uniform sampler2D u_gbuffer_depth;

uniform mat4 u_shadow_vp;
uniform sampler2D u_shadowmap;

uniform mat4 u_shadow_vp1;
uniform sampler2D u_shadowmap1;

uniform vec3 u_camera_position;

out vec4 FragColor;

// Lighing
uniform vec3 u_ambient_light;

const int MAX_LIGHT_COUNT = 10;
uniform int u_light_count;
uniform int u_light_type[MAX_LIGHT_COUNT];
uniform vec3 u_light_positions[MAX_LIGHT_COUNT];
uniform vec3 u_light_colors[MAX_LIGHT_COUNT];
uniform float u_light_intensities[MAX_LIGHT_COUNT];

float get_shadow_depth(vec3 world_pos) {
	vec4 fragment_shadow = u_shadow_vp * vec4(world_pos, 1.0);
	fragment_shadow = fragment_shadow / fragment_shadow.w;
	vec2 frag_shadows_uv = ((fragment_shadow.xy) * 0.5) + vec2(0.5);

	if (frag_shadows_uv.x < 0.0 && frag_shadows_uv.x > 1.0 && frag_shadows_uv.y < 0.0 && frag_shadows_uv.y > 1.0) {
		return 0.0;
	}

	float shadow_map_depth = texture(u_shadowmap, frag_shadows_uv).x;
	float frag_depth = ((((fragment_shadow.z - 0.0001))) * 0.5) + 0.5;
	if (shadow_map_depth < frag_depth) {
		return 0.0;
	}

	return 1.0;
}

float get_shadow_depth1(vec3 world_pos) {
	vec4 fragment_shadow = u_shadow_vp1 * vec4(world_pos, 1.0);
	fragment_shadow = fragment_shadow / fragment_shadow.w;
	vec2 frag_shadows_uv = ((fragment_shadow.xy) * 0.5) + vec2(0.5);

	if (frag_shadows_uv.x < 0.0 && frag_shadows_uv.x > 1.0 && frag_shadows_uv.y < 0.0 && frag_shadows_uv.y > 1.0) {
		return 0.0;
	}

	float shadow_map_depth = texture(u_shadowmap1, frag_shadows_uv).x;
	float frag_depth = ((((fragment_shadow.z - 0.0001))) * 0.5) + 0.5;
	if (shadow_map_depth < frag_depth) {
		return 0.0;
	}

	return 1.0;
}

void main()
{
	vec2 uv = v_uv;
	
	vec2 uv_clip = uv * 2.0 -1.0;

	//uv.y = 1.0 - uv.y;

	float depth = texture(u_gbuffer_depth, uv).r * 2.0 - 1.0;
	vec4 color = texture(u_gbuffer_albedo, uv);
	vec3 N = normalize(texture(u_gbuffer_normal, uv).rgb * 2.0 - 1.0);
	float normal_w = texture(u_gbuffer_normal, uv).w;

	vec4 not_norm_world_pos = u_inv_vp_mat * vec4(uv_clip.x, uv_clip.y, depth, 1.0);
	vec3 world_pos = not_norm_world_pos.xyz / not_norm_world_pos.w;

	// Ambient contributions
	vec3 outgoing_light = u_ambient_light;

	vec3 V = normalize(u_camera_position - world_pos);

	// Evaluate light contribution
	for(int i = 0; i < u_light_count; i++) {
		vec3 L = normalize(u_light_positions[i] - world_pos);
		vec3 R = reflect(-L, N);

		float light_dist = distance(u_light_positions[i], world_pos);
		vec3 light_attenuation = (u_light_intensities[i] * u_light_colors[i]) / (1.0+(light_dist*light_dist));

		float shadow = get_shadow_depth(world_pos);

		// Diffuse contribution
		if (i == 3) {
			outgoing_light += clamp(dot(L, N), 0.0, 1.0) * light_attenuation * shadow;
			outgoing_light += pow(clamp(dot(R, V), 0.0, 1.0), 64.0) * light_attenuation * shadow;
		} else {
			outgoing_light += clamp(dot(L, N), 0.0, 1.0) * light_attenuation;
			outgoing_light += pow(clamp(dot(R, V), 0.0, 1.0), 64.0) * light_attenuation;
		}
		
	}

	if (depth == 1.0) {
		outgoing_light = color.xyz;
	}

	// resulting_color = (ambeint + diffuse + specular) * base_color
	color *= vec4(outgoing_light, 1.0);

	FragColor =color;
}

\phong.fs

#version 330 core

in vec3 v_position;
in vec3 v_world_position;
in vec3 v_normal;
in vec2 v_uv;
in vec4 v_color;

uniform vec4 u_color;
uniform sampler2D u_texture;
uniform float u_time;
uniform float u_alpha_cutoff;
uniform vec3 u_camera_position;

uniform mat4 u_shadow_vp;
uniform sampler2D u_shadowmap;

uniform mat4 u_shadow_vp1;
uniform sampler2D u_shadowmap1;

out vec4 FragColor;

// Lighing
uniform vec3 u_ambient_light;

const int MAX_LIGHT_COUNT = 10;
uniform int u_light_count;
uniform int u_light_type[MAX_LIGHT_COUNT];
uniform vec3 u_light_positions[MAX_LIGHT_COUNT];
uniform vec3 u_light_dirs[MAX_LIGHT_COUNT];
uniform vec3 u_light_colors[MAX_LIGHT_COUNT];
uniform float u_light_intensities[MAX_LIGHT_COUNT];
uniform vec2 u_cone_data[MAX_LIGHT_COUNT];

float get_shadow_depth(vec3 world_pos) {
	vec4 fragment_shadow = u_shadow_vp * vec4(world_pos, 1.0);
	fragment_shadow = fragment_shadow / fragment_shadow.w;
	vec2 frag_shadows_uv = ((fragment_shadow.xy) * 0.5) + vec2(0.5);

	if (frag_shadows_uv.x < 0.0 && frag_shadows_uv.x > 1.0 && frag_shadows_uv.y < 0.0 && frag_shadows_uv.y > 1.0) {
		return 1.0;
	}

	float shadow_map_depth = texture(u_shadowmap, frag_shadows_uv).x;
	float frag_depth = ((((fragment_shadow.z - 0.0001))) * 0.5) + 0.5;
	if (shadow_map_depth < frag_depth) {
		return 1.0;
	}

	return 1.0;
}

float get_shadow_depth1(vec3 world_pos) {
	vec4 fragment_shadow = u_shadow_vp1 * vec4(world_pos, 1.0);
	fragment_shadow = fragment_shadow / fragment_shadow.w;
	vec2 frag_shadows_uv = (vec2(fragment_shadow.x, fragment_shadow.y) * 0.5) + vec2(0.5);

	if (frag_shadows_uv.x < 0.0 && frag_shadows_uv.x > 1.0 && frag_shadows_uv.y < 0.0 && frag_shadows_uv.y > 1.0) {
		return 1.0;
	}

	float shadow_map_depth = texture(u_shadowmap1, frag_shadows_uv).x;
	float frag_depth = ((((fragment_shadow.z - 0.0001))) * 0.5) + 0.5;
	if (shadow_map_depth < frag_depth) {
		return 1.0;
	}

	return 1.0;
}

void main()
{
	vec2 uv = v_uv;
	// Base color
	vec4 color = u_color;
	color *= texture2D( u_texture, uv );

	if(color.a < u_alpha_cutoff)
		discard;

	vec3 N = normalize(v_normal);

	// Ambient contributions
	vec3 outgoing_light = u_ambient_light;

	vec3 V = normalize(u_camera_position - v_world_position);

	// Evaluate light contribution
	for(int i = 0; i < 10; i++) {
		if (u_light_count <= i) {
			break;
		}
		vec3 L = normalize(u_light_positions[i] - v_world_position);

		float light_dist = distance(u_light_positions[i], v_world_position);
		vec3 light_attenuation = (u_light_intensities[i] * u_light_colors[i]) / (1.0+(light_dist*light_dist));;

		if (u_light_type[i] == 3) {
			L = normalize(u_light_dirs[i]);
			light_attenuation = u_light_intensities[i] * u_light_colors[i] * get_shadow_depth(v_world_position);
		} else if (u_light_type[i] == 2) {
			vec2 cone_data = u_cone_data[i];
			float minus_l_dot_d = clamp(dot(L, normalize(u_light_dirs[i])), 0.0, 1.0);

			if (minus_l_dot_d >= (cone_data.x)) {
				//light_attenuation *= vec3(pow(minus_l_dot_d, cone_data.y * 200.0));
				light_attenuation *= clamp((minus_l_dot_d - cos(cone_data.y)) / (cos(cone_data.x)- cos(cone_data.y)), 0.0, 1.0) * get_shadow_depth1(v_world_position);

				
			} else {
				light_attenuation = vec3(0.0);
			}

			//light_attenuation = vec3(minus_l_dot_d);
		}

		vec3 R = reflect(-L, N);


		// Diffuse contribution
		if (i == 3) {
			outgoing_light += clamp(dot(L, N), 0.0, 1.0) * light_attenuation;
			outgoing_light += pow(clamp(dot(R, V), 0.0, 1.0), 64.0) * light_attenuation;
		} else {
			outgoing_light += clamp(dot(L, N), 0.0, 1.0) * light_attenuation;
			outgoing_light += pow(clamp(dot(R, V), 0.0, 1.0), 64.0) * light_attenuation;
		}
		
	}

	// resulting_color = (ambeint + diffuse + specular) * base_color
	color *= vec4(outgoing_light, 1.0);

	FragColor = color;
}

\skybox.fs

#version 330 core

in vec3 v_position;
in vec3 v_world_position;

uniform samplerCube u_texture;
uniform vec3 u_camera_position;

layout(location = 0) out vec4 gbuffer_albedo;
layout(location = 1) out vec4 gbuffer_normal_mat;
void main()
{
	vec3 E = v_world_position - u_camera_position;
	vec4 color = texture( u_texture, E );
	gbuffer_albedo = color;
	gbuffer_normal_mat = vec4(0.0);
}


\multi.fs

#version 330 core

in vec3 v_position;
in vec3 v_world_position;
in vec3 v_normal;
in vec2 v_uv;

uniform vec4 u_color;
uniform sampler2D u_texture;
uniform float u_time;
uniform float u_alpha_cutoff;

layout(location = 0) out vec4 FragColor;
layout(location = 1) out vec4 NormalColor;

void main()
{
	vec2 uv = v_uv;
	vec4 color = u_color;
	color *= texture( u_texture, uv );

	if(color.a < u_alpha_cutoff)
		discard;

	vec3 N = normalize(v_normal);

	FragColor = color;
	NormalColor = vec4(N,1.0);
}


\depth.fs

#version 330 core

uniform vec2 u_camera_nearfar;
uniform sampler2D u_texture; //depth map
in vec2 v_uv;
out vec4 FragColor;

void main()
{
	float n = u_camera_nearfar.x;
	float f = u_camera_nearfar.y;
	float z = texture2D(u_texture,v_uv).x;
	if( n == 0.0 && f == 1.0 )
		FragColor = vec4(z);
	else
		FragColor = vec4( n * (z + 1.0) / (f + n - z * (f - n)) );
}


\instanced.vs

#version 330 core

in vec3 a_vertex;
in vec3 a_normal;
in vec2 a_coord;

in mat4 u_model;

uniform vec3 u_camera_pos;

uniform mat4 u_viewprojection;

//this will store the color for the pixel shader
out vec3 v_position;
out vec3 v_world_position;
out vec3 v_normal;
out vec2 v_uv;

void main()
{	
	//calcule the normal in camera space (the NormalMatrix is like ViewMatrix but without traslation)
	v_normal = (u_model * vec4( a_normal, 0.0) ).xyz;
	
	//calcule the vertex in object space
	v_position = a_vertex;
	v_world_position = (u_model * vec4( a_vertex, 1.0) ).xyz;
	
	//store the texture coordinates
	v_uv = a_coord;

	//calcule the position of the vertex using the matrices
	gl_Position = u_viewprojection * vec4( v_world_position, 1.0 );
}


\light_volume.fs
#version 330 core

in vec3 v_position;
in vec3 v_world_position;
in vec3 v_normal;
in vec2 v_uv;
in vec4 v_color;

uniform mat4 u_inv_vp_mat;

uniform sampler2D u_gbuffer_albedo;
uniform sampler2D u_gbuffer_normal;
uniform sampler2D u_gbuffer_depth;

uniform mat4 u_shadow_vp;
uniform sampler2D u_shadowmap;

uniform mat4 u_shadow_vp1;
uniform sampler2D u_shadowmap1;

uniform vec3 u_camera_position;

uniform vec2 u_res_inv;

out vec4 FragColor;

#include "PBR"

// Lighing
uniform vec3 u_ambient_light;

const int MAX_LIGHT_COUNT = 10;
uniform int u_light_count;
uniform int u_light_type[MAX_LIGHT_COUNT];
uniform vec3 u_light_positions[MAX_LIGHT_COUNT];
uniform vec3 u_light_colors[MAX_LIGHT_COUNT];
uniform float u_light_intensities[MAX_LIGHT_COUNT];
uniform vec2 u_cone_data[MAX_LIGHT_COUNT];
uniform vec3 u_light_dirs[MAX_LIGHT_COUNT];


uniform int u_index;

float get_shadow_depth(vec3 world_pos) {
	vec4 fragment_shadow = u_shadow_vp * vec4(world_pos, 1.0);
	fragment_shadow = fragment_shadow / fragment_shadow.w;
	vec2 frag_shadows_uv = ((fragment_shadow.xy) * 0.5) + vec2(0.5);

	if (frag_shadows_uv.x < 0.0 && frag_shadows_uv.x > 1.0 && frag_shadows_uv.y < 0.0 && frag_shadows_uv.y > 1.0) {
		return 0.0;
	}

	float shadow_map_depth = texture(u_shadowmap, frag_shadows_uv).x;
	float frag_depth = ((((fragment_shadow.z - 0.0001))) * 0.5) + 0.5;
	if (shadow_map_depth < frag_depth) {
		return 0.0;
	}

	return 1.0;
}

void main()
{
	vec2 uv = gl_FragCoord.xy * u_res_inv;
	
	float depth = texture(u_gbuffer_depth, uv).r * 2.0 - 1.0;
	vec4 color = texture(u_gbuffer_albedo, uv);
	vec3 N = normalize(texture(u_gbuffer_normal, uv).rgb * 2.0 - 1.0);
	float normal_w = texture(u_gbuffer_normal, uv).w;

	vec3 albedo = color.xyz;
	float roughness = color.w;
	float metalness = normal_w;

	vec2 uv_clip = uv * 2.0 -1.0;

	vec4 not_norm_world_pos = u_inv_vp_mat * vec4(uv_clip.x, uv_clip.y, depth, 1.0);
	vec3 world_pos = not_norm_world_pos.xyz / not_norm_world_pos.w;

	// Ambient contributions
	vec3 outgoing_light = vec3(0.0);

	vec3 V = normalize(u_camera_position - world_pos);

	// Evaluate light contribution
		vec3 L = normalize(u_light_positions[u_index] - world_pos);

		float light_dist = distance(u_light_positions[u_index], world_pos);
		vec3 light_attenuation = (u_light_intensities[u_index] * u_light_colors[u_index]) / (1.0+(light_dist*light_dist));

		if (u_light_type[u_index] == 2) {
			vec2 cone_data = u_cone_data[u_index];
			float minus_l_dot_d = clamp(dot(L, normalize(u_light_dirs[u_index])), 0.0, 1.0);

			if (minus_l_dot_d >= (cone_data.x)) {
				//light_attenuation *= vec3(pow(minus_l_dot_d, cone_data.y * 200.0));
				light_attenuation *= clamp((minus_l_dot_d - cos(cone_data.y)) / (cos(cone_data.x)- cos(cone_data.y)), 0.0, 1.0);

				
			} else {
				light_attenuation = vec3(0.0);
			}

			//light_attenuation = vec3(minus_l_dot_d);
		}

		vec3 R = reflect(-L, N);
		vec3 H = normalize(V + L);

		vec3 direct_light = brdf_cook_torrance(albedo, roughness, metalness, N, H, L, V);
		// Diffuse contribution
		outgoing_light += clamp(dot(L, N), 0.0, 1.0) * light_attenuation * direct_light;

	// resulting_color = (ambeint + diffuse + specular) * base_color
	color *= vec4(outgoing_light, 1.0);

	FragColor = color;
}

\ssao_pass.fs
#version 330 core

in vec2 v_uv;

uniform mat4 u_v_mat;
uniform mat4 u_p_mat;
uniform mat4 u_inv_p_mat;

uniform sampler2D u_gbuffer_normal;
uniform sampler2D u_gbuffer_depth;

uniform vec2 u_res_inv;

const int MAX_SAMPLE_COUNT = 30;
uniform int u_sample_count;
uniform vec3 u_sample_pos[MAX_SAMPLE_COUNT];
uniform float u_sample_radius;

out vec4 FragColor;

vec3 points[16] = {
vec3(0.355647266, 0.0937004983, -0.0243267361),
vec3(-0.340726703, -2.97872749e-08, -0.300947130),
vec3(-0.469363153, 0.0652965903, -0.140362382),
vec3(-0.272493631, 0.0627098009, 0.125692934),
vec3(-0.0723427683, 0.366701156, 0.323128462),
vec3(-0.0691736117, 0.00391585613, 0.486232638),
vec3(0.177007183, -0.336863309, -0.0491127186),
vec3(-0.146231636, -0.111964323, 0.249020889),
vec3(-0.0991556719, 0.366977900, -0.0923202634),
vec3(-0.0922925398, -0.425431788, 0.210904568),
vec3(-0.390294969, 0.134999231, 0.209265247),
vec3(-0.318935096, -0.0776300579, 0.358989298),
vec3(-0.115982905, -0.201863468, -0.102212638),
vec3(0.337634087, -0.265310466, 0.190767586),
vec3(-0.183409393, 0.377381384, -0.148681372),
vec3(-0.161931634, -0.269925982, -0.324702382)
};


vec3 get_view_pos_from_UVs(vec2 uv, float z_delta) {
	float depth = texture(u_gbuffer_depth, uv).r;

	vec4 clip_coords = vec4(uv.x * 2.0 - 1.0, uv.y * 2.0 - 1.0, depth* 2.0 - 1.0, 1.0);

	vec4 not_norm_view_coord = u_inv_p_mat * clip_coords;

	return (not_norm_view_coord.xyz + vec3(0.0, 0.0, z_delta)) / not_norm_view_coord.w;
}

void main()
{
	vec2 uv = gl_FragCoord.xy * u_res_inv + 0.5 * u_res_inv;

	float sample_radius = 0.08;
	
	float depth = texture(u_gbuffer_depth, uv).r;

	if (depth >= 1.0) {
		FragColor = vec4(1.0);
		return;
	}

	vec3 view_pos = get_view_pos_from_UVs(uv, 0.0);

	float ao_term = 0.0;
	float sample_count = 0.0;

	for(int i = 0; i < 30; i++) {
		vec3 sample_pos = u_sample_pos[i] * sample_radius + view_pos;

		//sample_pos += normal_view  * 0.0001;

		vec4 proj_sample = u_p_mat * vec4(sample_pos, 1.0);
		proj_sample /= proj_sample.w;

		vec2 uv_clip = proj_sample.xy * 0.5 + 0.5;

		float depth = texture(u_gbuffer_depth, uv_clip).r;

		float clip_depth = depth * 2.0 - 1.0;

		vec4 pr = u_inv_p_mat * vec4(uv_clip, clip_depth, 1.0);
		pr /= pr.w;

		sample_count += 1.0;


		if (clip_depth > (proj_sample.z)) {
			ao_term += 1.0;// * smoothstep(0.0, 1.0, sample_radius / abs(sample_pos.z - pr.z));
		}
	}

	ao_term /= sample_count;

	FragColor = vec4(ao_term);
}


\first_pass.fs
#version 330 core

in vec2 v_uv;

uniform mat4 u_inv_vp_mat;

uniform sampler2D u_gbuffer_albedo;
uniform sampler2D u_gbuffer_normal;
uniform sampler2D u_gbuffer_depth;


uniform mat4 u_shadow_vp;
uniform sampler2D u_shadowmap;

uniform mat4 u_shadow_vp1;
uniform sampler2D u_shadowmap1;

uniform vec3 u_camera_position;

uniform vec2 u_res_inv;

out vec4 FragColor;

// Lighing
uniform vec3 u_ambient_light;

uniform vec3 u_light_dir;
uniform vec3 u_light_color;
uniform float u_light_intensity;

uniform int u_index;

#include "PBR"

float get_shadow_depth(vec3 world_pos) {
	vec4 fragment_shadow = u_shadow_vp * vec4(world_pos, 1.0);
	fragment_shadow = fragment_shadow / fragment_shadow.w;
	vec2 frag_shadows_uv = ((fragment_shadow.xy) * 0.5) + vec2(0.5);

	if (frag_shadows_uv.x < 0.0 && frag_shadows_uv.x > 1.0 && frag_shadows_uv.y < 0.0 && frag_shadows_uv.y > 1.0) {
		return 0.0;
	}

	float shadow_map_depth = texture(u_shadowmap, frag_shadows_uv).x;
	float frag_depth = ((((fragment_shadow.z - 0.001))) * 0.5) + 0.5;
	if (shadow_map_depth < frag_depth) {
		return 0.0;
	}

	return 1.0;
}

void main()
{
	vec2 uv = gl_FragCoord.xy * u_res_inv;
	
	float depth = texture(u_gbuffer_depth, uv).r * 2.0 - 1.0;
	vec4 color = texture(u_gbuffer_albedo, uv);
	vec3 N = normalize(texture(u_gbuffer_normal, uv).rgb * 2.0 - 1.0);
	float normal_w = texture(u_gbuffer_normal, uv).w;

	vec3 albedo = color.xyz;
	float roughness = color.w;
	float metalness = normal_w;

	vec2 uv_clip = uv * 2.0 -1.0;

	vec4 not_norm_world_pos = u_inv_vp_mat * vec4(uv_clip.x, uv_clip.y, depth, 1.0);
	vec3 world_pos = not_norm_world_pos.xyz / not_norm_world_pos.w;

	// Ambient contributions
	vec3 outgoing_light = u_ambient_light;

	vec3 V = normalize(u_camera_position - world_pos);

	// Evaluate light contribution
		vec3 L = normalize(u_light_dir);
		vec3 H = normalize(V + L);

		vec3 light_attenuation = (u_light_intensity * u_light_color);

		float shadow = get_shadow_depth(world_pos);

		vec3 direct_light = brdf_cook_torrance(albedo, roughness, metalness, N, H, L, V);
		outgoing_light += clamp(dot(L, N), 0.0, 1.0) * direct_light * light_attenuation * shadow;		

	// resulting_color = (ambeint + diffuse + specular) * base_color
	color *= vec4(outgoing_light, 1.0);

	FragColor = color;
}