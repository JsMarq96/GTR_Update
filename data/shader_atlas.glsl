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

phong basic.vs phong.fs

\PBR

vec3 F_Schlick(vec3 v, vec3 h, vec3 f0) {
	return f0 + (1.0 - f0) * pow(1.0 - dot(h, v), 5.0);
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

	if(color.a < u_alpha_cutoff)
		discard;

	vec3 N = normalize(v_normal);

	gbuffer_albedo = color;
	gbuffer_normal_mat = vec4(N * 0.5 + 0.5, 64.0);
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


		// Diffuse contribution
		outgoing_light += clamp(dot(L, N), 0.0, 1.0) * light_attenuation;
		outgoing_light += pow(clamp(dot(R, V), 0.0, 1.0), 64.0) * light_attenuation;

	// resulting_color = (ambeint + diffuse + specular) * base_color
	color *= vec4(outgoing_light, 1.0);

	FragColor = color;
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

	vec2 uv_clip = uv * 2.0 -1.0;

	vec4 not_norm_world_pos = u_inv_vp_mat * vec4(uv_clip.x, uv_clip.y, depth, 1.0);
	vec3 world_pos = not_norm_world_pos.xyz / not_norm_world_pos.w;

	// Ambient contributions
	vec3 outgoing_light = u_ambient_light;

	vec3 V = normalize(u_camera_position - world_pos);

	// Evaluate light contribution
		vec3 L = normalize(u_light_dir);
		vec3 R = reflect(-L, N);

		vec3 light_attenuation = (u_light_intensity * u_light_color);

		float shadow = get_shadow_depth(world_pos);

		// Diffuse contribution
		outgoing_light += clamp(dot(L, N), 0.0, 1.0) * light_attenuation * shadow;
		outgoing_light += pow(clamp(dot(R, V), 0.0, 1.0), 64.0) * light_attenuation * shadow;
		

	// resulting_color = (ambeint + diffuse + specular) * base_color
	color *= vec4(outgoing_light, 1.0);

	FragColor = color;
}