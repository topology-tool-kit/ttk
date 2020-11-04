#pragma once

// VTK Module
#include <ttkCinemaDarkroomModule.h>
#include <ttkAlgorithm.h>
#include <vtkSmartPointer.h>

#include <unordered_map>

class vtkImageData;
class vtkPolyData;
class vtkActor;
class vtkRenderer;
class vtkRenderWindow;
class vtkWindowToImageFilter;

class TTKCINEMADARKROOM_EXPORT ttkCinemaDarkroomShader : public ttkAlgorithm {

private:
  struct Replacement {
    std::vector<double> values;
    bool isInt{false};

    Replacement(const std::vector<double>& values_, const bool& isInt_)
      : values(values_), isInt(isInt_)
    {
    }

    std::string toString() const {
      std::string result = "";
      if(this->values.size() == 0){
        return "";
      }

      if(this->values.size()>1){
        if(this->isInt)
          result += "i";

        result += "vec"+std::to_string(this->values.size())+"(";
      }

      if(this->isInt)
        result += std::to_string((int)this->values[0]);
      else
        result += std::to_string(this->values[0]);

      for(size_t i=1; i<this->values.size(); i++)
        if(this->isInt)
          result += "," + std::to_string((int)this->values[i]);
        else
          result += "," + std::to_string(this->values[i]);

      if(this->values.size()>1)
        result += ")";

      return result;
    }
  };

  std::unordered_map<std::string, Replacement> Replacements;

  // Quad
  vtkSmartPointer<vtkPolyData> FullScreenQuad;
  vtkSmartPointer<vtkActor> FullScreenQuadActor;

  // Rendering Tools
  vtkSmartPointer<vtkRenderer> Renderer;
  vtkSmartPointer<vtkRenderWindow> RenderWindow;
  vtkSmartPointer<vtkWindowToImageFilter> RenderWindowToImageFilter;

public:

  static ttkCinemaDarkroomShader *New();
  vtkTypeMacro(ttkCinemaDarkroomShader, ttkAlgorithm);

protected:
  ttkCinemaDarkroomShader();
  ~ttkCinemaDarkroomShader() override;

  std::string PerformReplacements(const std::string& source);
  int AddReplacement(const std::string& name, const std::vector<double>& values, const bool& isInt = false);

  int CreateFullScreenQuad();
  int CreateRenderer();
  int InitRenderer(vtkImageData* outputImage);
  int SetTexture(vtkImageData* image, int arrayIdx, int textureIdx, int* textureProperties);

  virtual std::string GetVertexShaderCode();
  virtual std::string GetFragmentShaderCode();
  virtual int Render(vtkImageData* image, const std::string& name);

  int FillInputPortInformation(int port, vtkInformation *info) override;
  int FillOutputPortInformation(int port, vtkInformation *info) override;
};