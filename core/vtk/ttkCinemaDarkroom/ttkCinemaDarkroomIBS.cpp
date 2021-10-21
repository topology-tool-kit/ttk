#include <ttkCinemaDarkroomIBS.h>

#include <vtkImageData.h>
#include <vtkInformation.h>
#include <vtkObjectFactory.h>

vtkStandardNewMacro(ttkCinemaDarkroomIBS);

ttkCinemaDarkroomIBS::ttkCinemaDarkroomIBS() : ttkCinemaDarkroomShader() {
  this->setDebugMsgPrefix("CinemaDarkroomIBS");
}

ttkCinemaDarkroomIBS::~ttkCinemaDarkroomIBS() {
}

std::string ttkCinemaDarkroomIBS::GetFragmentShaderCode() {
  return std::string(R"(
//VTK::System::Dec // always start with these lines in your FS
//VTK::Output::Dec // always start with these lines in your FS

// #extension GL_OES_standard_derivatives : enable

varying vec4 vPos;

uniform sampler2D tex0; // color
uniform sampler2D tex1; // depth
uniform sampler2D tex2; // ao

float readDepth( const in vec2 coord ){
    return texture2D( tex1, coord ).r;
}

void main() {
    vec3 color = texture2D( tex0, vPos.xy ).rgb;
    float ao = texture2D( tex2, vPos.xy ).r;
    float depth = readDepth(vPos.xy);

    // Compute Luminance
    vec3 lumcoeff = vec3( 0.299, 0.587, 0.114 );
    vec3 luminance = vec3( dot( color, lumcoeff ) );

    // Silhouette Effect
    vec2 pixelSize = 1./cResolution;
    vec3 eps = 2.0*vec3( pixelSize.x, pixelSize.y, 0 );
    float depthN = readDepth(vPos.xy + eps.zy);
    float depthE = readDepth(vPos.xy + eps.xz);
    float depthS = readDepth(vPos.xy - eps.zy);
    float depthW = readDepth(vPos.xy - eps.xz);

    float dxdz = abs(depthE-depthW);
    float dydz = abs(depthN-depthS);
    // float dxdz = dFdx(depth);
    // float dydz = dFdy(depth);

    vec3 n = normalize( vec3(dxdz, dydz, 1./cStrength) );
    vec3 lightPos = vec3(0,0,1);
    float lightInt = 1.0*dot(n,normalize(lightPos));

    vec3 outputColor = vec3( color * mix( vec3(ao), vec3(1.0), luminance * cLuminance ) );

    outputColor = outputColor*cAmbient + outputColor*lightInt;

    gl_FragColor = vec4(outputColor, depth>0.99 ? 0.0 : 1.0);
}
  )");
}

int ttkCinemaDarkroomIBS::RequestData(vtkInformation *ttkNotUsed(request),
                                      vtkInformationVector **inputVector,
                                      vtkInformationVector *outputVector) {

  auto inputImage = vtkImageData::GetData(inputVector[0]);
  auto outputImage = vtkImageData::GetData(outputVector);
  outputImage->ShallowCopy(inputImage);

  this->InitRenderer(outputImage);

  this->AddReplacement("cStrength", {this->Strength});
  this->AddReplacement("cLuminance", {this->Luminance});
  this->AddReplacement("cAmbient", {this->Ambient});

  if(!this->AddTexture(outputImage, 0, 0))
    return 0;
  if(!this->AddTexture(outputImage, 1, 1))
    return 0;
  if(!this->AddTexture(outputImage, 2, 2))
    return 0;

  this->Render(outputImage, "IBS");

  return 1;
}