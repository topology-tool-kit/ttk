/// \ingroup vtk
/// \class ttkCinemaDarkroomShader
/// \author Jonas Lukasczyk <jl@jluk.de>
/// \date 01.11.2020
///
/// \brief Base Class for all CinemaDarkroom Shaders.
///
/// \param Input vtkImageData.
/// \param Output vtkImageData.
///
/// This class provides all utility functions to create and use an OpenGL render
/// context. The core concept of this class is to make a shallow copy of its
/// vtkImageData input, and then add a point data array to the output that
/// corresponds to a render pass. This class provides various abstract concepts
/// to configure this render pass in the RequestData method:
///
/// 1) Subclasses can override the GetVertexShaderCode and the
/// GetFragmentShaderCode methods to provide their own shaders. By default, the
/// vertex shader renders a fullscreen quad, and the fragment shader simply
/// produces a red color.
///
/// 2) Instead of actually passing shader uniforms, this class uses string
/// replacements to feed shader parameters into a render pass. One can add a
/// replacement with the AddReplacement method, where type conversion to GLSL is
/// done automatically.
///
/// 3) Point data arrays of the input can be passed into a render pass as data
/// textures via the AddTexture method.
///
/// \b Related \b Publication:
/// "Cinema Darkroom: A Deferred Rendering Framework for Large-Scale Datasets".
/// J. Lukasczyk, C. Garth, M. Larsen, W. Engelke, I. Hotz, D. Rogers, J.
/// Ahrens, and R. Maciejewski. IEEE 10th Symposium on Large Data Analysis and
/// Visualization (LDAV), 2020.
///
/// \sa ttkAlgorithm

#pragma once

// VTK Module
#include <ttkAlgorithm.h>
#include <ttkCinemaDarkroomModule.h>
#include <vtkSmartPointer.h>

#include <unordered_map>

class vtkImageData;
class vtkPolyData;
class vtkActor;
class vtkRenderer;
class vtkRenderWindow;

class TTKCINEMADARKROOM_EXPORT ttkCinemaDarkroomShader : public ttkAlgorithm {

private:
  struct Replacement {
    std::vector<double> values;
    bool isInt{false};

    Replacement(const std::vector<double> &values_, const bool &isInt_)
      : values(values_), isInt(isInt_) {
    }

    std::string toString() const {
      std::string result = "";
      if(this->values.size() == 0) {
        return "";
      }

      if(this->values.size() > 1) {
        if(this->isInt)
          result += "i";

        result += "vec" + std::to_string(this->values.size()) + "(";
      }

      if(this->isInt)
        result += std::to_string((int)this->values[0]);
      else
        result += std::to_string(this->values[0]);

      for(size_t i = 1; i < this->values.size(); i++)
        if(this->isInt)
          result += "," + std::to_string((int)this->values[i]);
        else
          result += "," + std::to_string(this->values[i]);

      if(this->values.size() > 1)
        result += ")";

      return result;
    }
  };

  std::unordered_map<std::string, Replacement> Replacements;

  vtkSmartPointer<vtkPolyData> FullScreenQuad;
  vtkSmartPointer<vtkActor> FullScreenQuadActor;
  vtkSmartPointer<vtkRenderer> Renderer;
  vtkSmartPointer<vtkRenderWindow> RenderWindow;

public:
  static ttkCinemaDarkroomShader *New();
  vtkTypeMacro(ttkCinemaDarkroomShader, ttkAlgorithm);

protected:
  ttkCinemaDarkroomShader();
  ~ttkCinemaDarkroomShader() override;

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;

  /// Replaces all palceholder strings of the input string based on previously
  /// added replacements.
  std::string PerformReplacements(const std::string &input);

  /// Adds or updates a replacement that will be later used during the render
  /// pass to update shader strings.
  int AddReplacement(const std::string &name,
                     const std::vector<double> &values,
                     const bool &isInt = false);

  int CreateFullScreenQuad();
  int CreateRenderer();

  /// Initializes the render context and updates the internal render-window size
  /// based on the dimension of the output image.
  int InitRenderer(vtkImageData *outputImage);

  /// Adds a texture to the render pass with the name tex{textureIdx} whose
  /// values correspond to the point data array specified via
  /// SetInputArrayToProcess(arrayIdx).
  int AddTexture(vtkImageData *image, int arrayIdx, int textureIdx);

  virtual std::string GetVertexShaderCode();
  virtual std::string GetFragmentShaderCode();

  /// Performs a single render pass and adds the result to the image as a point
  /// data array with the specified name.
  virtual int Render(vtkImageData *image, const std::string &name);
};