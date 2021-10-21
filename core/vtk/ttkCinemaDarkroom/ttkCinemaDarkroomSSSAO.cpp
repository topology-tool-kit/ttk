#include <ttkCinemaDarkroomSSSAO.h>

#include <vtkImageData.h>
#include <vtkInformation.h>
#include <vtkObjectFactory.h>

vtkStandardNewMacro(ttkCinemaDarkroomSSSAO);

ttkCinemaDarkroomSSSAO::ttkCinemaDarkroomSSSAO() : ttkCinemaDarkroomShader() {
  this->setDebugMsgPrefix("CinemaDarkroomSSSAO");
}

ttkCinemaDarkroomSSSAO::~ttkCinemaDarkroomSSSAO() {
}

std::string ttkCinemaDarkroomSSSAO::GetFragmentShaderCode() {
  return std::string(R"(
//VTK::System::Dec // always start with these lines in your FS
//VTK::Output::Dec // always start with these lines in your FS

// -----------------------------------------------------------------------------
// Custom SSAO implementation based on Blender's Viewport SSAO (https://github.com/blender/blender)
// Roots can be traced back to Arkano22 (https://www.gamedev.net/forums/topic/550699-ssao-no-halo-artifacts/)
// and Martins Upitis (martinsh) (http://devlog-martinsh.blogspot.com/search/label/SSAO)
// -----------------------------------------------------------------------------

uniform sampler2D tex0;

varying vec4 vPos;

#define DL 2.399963229728653  // PI * ( 3.0 - sqrt( 5.0 ) )
#define EULER 2.718281828459045

float readDepth( const in vec2 coord ){
    return texture2D( tex0, coord ).r;
}

const float gDisplace = 0.5;  // gauss bell center
float compareDepths( const in float depth1, const in float depth2, inout int far ) {
    float garea = 16.0;        // gauss bell width
    float diff = ( depth1 - depth2 ) * 100.0;  // depth difference (0-100)

    // reduce left bell width to avoid self-shadowing
    if(diff<gDisplace){
        garea = cDiffArea;
    } else {
        far = 1;
    }

    float dd = diff - gDisplace;
    return pow( EULER, -2.0 * ( dd * dd ) / ( garea * garea ) );
}

float calcAO( float depth, float dw, float dh, vec2 uv ) {
    vec2 vv = vec2( dw, dh );
    vec2 coord1 = uv + vv;
    vec2 coord2 = uv - vv;
    float temp1 = 0.0;
    float temp2 = 0.0;
    int far = 0;

    temp1 = compareDepths( depth, readDepth( coord1 ), far );
    if ( far > 0 ) {
        temp2 = compareDepths( readDepth( coord2 ), depth, far );
        temp1 += ( 1.0 - temp1 ) * temp2;
    }
    return temp1;
}

void main() {
    float depth = readDepth( vPos.xy );

    const float samplesF = cSamples;
    float occlusion = 0.0;

    float dz = 1.0 / samplesF;
    float l = 0.0;
    float z = 1.0 - dz / 2.0;

    float aspect = cResolution.y/cResolution.x;

    for(int i=0; i<cSamples; i++){
        float r = sqrt( 1.0 - z ) * cRadius;
        float pw = cos( l ) * r;
        float ph = sin( l ) * r;
        occlusion += calcAO( depth, pw * aspect, ph, vPos.xy );
        z = z - dz;
        l = l + DL;
    }

    float ao = 1.-occlusion/samplesF;
    gl_FragColor = vec4(ao,ao,ao, 1 );
}
  )");
}

int ttkCinemaDarkroomSSSAO::RequestData(vtkInformation *ttkNotUsed(request),
                                        vtkInformationVector **inputVector,
                                        vtkInformationVector *outputVector) {

  auto inputImage = vtkImageData::GetData(inputVector[0]);
  auto outputImage = vtkImageData::GetData(outputVector);
  outputImage->ShallowCopy(inputImage);

  int dim[3];
  outputImage->GetDimensions(dim);

  this->InitRenderer(outputImage);

  this->AddReplacement("cSamples", {(double)this->Samples}, true);
  this->AddReplacement("cRadius", {this->Radius});
  this->AddReplacement("cDiffArea", {this->DiffArea});

  if(!this->AddTexture(outputImage, 0, 0))
    return 0;

  this->Render(outputImage, "SSSAO");

  return 1;
}