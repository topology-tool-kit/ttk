#include <ttkCinemaDarkroomSSAO.h>

#include <vtkImageData.h>
#include <vtkInformation.h>
#include <vtkObjectFactory.h>

vtkStandardNewMacro(ttkCinemaDarkroomSSAO);

ttkCinemaDarkroomSSAO::ttkCinemaDarkroomSSAO() : ttkCinemaDarkroomShader() {
  this->setDebugMsgPrefix("CinemaDarkroomSSAO");
}

ttkCinemaDarkroomSSAO::~ttkCinemaDarkroomSSAO() {
}

std::string ttkCinemaDarkroomSSAO::GetFragmentShaderCode() {
  return std::string(R"(
//VTK::System::Dec // always start with these lines in your FS
//VTK::Output::Dec // always start with these lines in your FS

uniform sampler2D tex0;
varying vec4 vPos;

float readDepth( const in vec2 coord ){
    return texture2D( tex0, coord ).r;
}

vec3 computePos( const in vec2 coord ){
    return vec3(
        coord,
        readDepth(coord)
    );
}

vec3 computeNormal(){
    vec2 pixelSize = 1./cResolution;

    vec3 eps = 2.0*vec3( pixelSize.x, pixelSize.y, 0 );
    float depthN = readDepth(vPos.xy + eps.zy);
    float depthE = readDepth(vPos.xy + eps.xz);
    float depthS = readDepth(vPos.xy - eps.zy);
    float depthW = readDepth(vPos.xy - eps.xz);

    float dzdx = (depthE - depthW) / 2.0;
    float dzdy = (depthN - depthS) / 2.0;

    return  normalize(vec3(-dzdx, -dzdy, 1.0));
}

const vec2 poissonDisk[32] = vec2[](
  vec2( -0.94201624,  -0.39906216 ),
  vec2(  0.94558609,  -0.76890725 ),
  vec2( -0.094184101, -0.92938870 ),
  vec2(  0.34495938,   0.29387760 ),
  vec2( -0.91588581,   0.45771432 ),
  vec2( -0.81544232,  -0.87912464 ),
  vec2( -0.38277543,   0.27676845 ),
  vec2(  0.97484398,   0.75648379 ),
  vec2(  0.44323325,  -0.97511554 ),
  vec2(  0.53742981,  -0.47373420 ),
  vec2( -0.26496911,  -0.41893023 ),
  vec2(  0.79197514,   0.19090188 ),
  vec2( -0.24188840,   0.99706507 ),
  vec2( -0.81409955,   0.91437590 ),
  vec2(  0.19984126,   0.78641367 ),
  vec2(  0.14383161,  -0.14100790 ),
  vec2( -0.44201624,  -0.29906216 ),
  vec2(  0.94558609,  -0.46890725 ),
  vec2( -0.194184101, -0.42938870 ),
  vec2(  0.24495938,   0.99387760 ),
  vec2( -0.31588581,   0.45771432 ),
  vec2( -0.81544232,  -0.87912464 ),
  vec2( -0.08277543,   0.87676845 ),
  vec2(  0.57484398,   0.55648379 ),
  vec2(  0.74323325,  -0.27511554 ),
  vec2(  0.44298431,  -0.47373420 ),
  vec2( -0.21196911,  -0.22893023 ),
  vec2(  0.79197514,   0.12020188 ),
  vec2( -0.11184840,   0.99706507 ),
  vec2( -0.4309955,   0.111437590 ),
  vec2(  0.12344126,   0.78641367 ),
  vec2(  0.2183161,   -0.89100790 )
);

void main(){
    vec3 pos = computePos(vPos.xy);
    vec3 n = computeNormal();

    float occlusion = 0.0;

    vec2 aspect = vec2(cResolution.y/cResolution.x, 1) * cRadius;

    for (int i = 0; i < 32; ++i){
        // get sample
        vec2 sampleTexCoord = vPos.xy + poissonDisk[i] * aspect;
        vec3 samplePos      = computePos(sampleTexCoord);
        vec3 sampleDir      = normalize(pos - samplePos);

        // distance between SURFACE-POSITION and SAMPLE-POSITION
        float VPdistSP = distance(pos, samplePos);

        // angle between SURFACE-NORMAL and SAMPLE-DIRECTION (vector from SURFACE-POSITION to SAMPLE-POSITION)
        float dotNS = max(dot(sampleDir, n), 0.0);

        // occlusion factor
        float a = 1.0 - smoothstep(cDiffArea, cDiffArea * 2.0, VPdistSP);

        // aggregate
        occlusion += a*dotNS;
    }

    float ao = 1.-occlusion/32.0;
    gl_FragData[0] = vec4( ao,ao,ao, 1 );
}
  )");
}

int ttkCinemaDarkroomSSAO::RequestData(vtkInformation *ttkNotUsed(request),
                                       vtkInformationVector **inputVector,
                                       vtkInformationVector *outputVector) {

  auto inputImage = vtkImageData::GetData(inputVector[0]);
  auto outputImage = vtkImageData::GetData(outputVector);
  outputImage->ShallowCopy(inputImage);

  this->InitRenderer(outputImage);

  this->AddReplacement("cRadius", {this->Radius});
  this->AddReplacement("cDiffArea", {this->DiffArea});

  if(!this->AddTexture(outputImage, 0, 0))
    return 0;

  this->Render(outputImage, "SSAO");

  return 1;
}