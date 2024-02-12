#include "rasterizer.h"

using namespace std;

namespace CGL {

  RasterizerImp::RasterizerImp(PixelSampleMethod psm, LevelSampleMethod lsm,
    size_t width, size_t height,
    unsigned int sample_rate) {
    this->psm = psm;
    this->lsm = lsm;
    this->width = width;
    this->height = height;
    this->sample_rate = sample_rate;

    sample_buffer.resize(width * height * sample_rate, Color::White);
  }

  // Used by rasterize_point and rasterize_line
  void RasterizerImp::fill_pixel(size_t x, size_t y, Color c) {
    // TODO: Task 2: You might need to this function to fix points and lines (such as the black rectangle border in test4.svg)
    // NOTE: You are not required to implement proper supersampling for points and lines
    // It is sufficient to use the same color for all supersamples of a pixel for points and lines (not triangles)

    for (int i = 0; i < sample_rate; i++) {
      sample_buffer[sample_rate * (y * width + x) + i] = c;
    }
  }

  // Rasterize a point: simple example to help you start familiarizing
  // yourself with the starter code.
  //
  void RasterizerImp::rasterize_point(float x, float y, Color color) {
    // fill in the nearest pixel
    int sx = (int)floor(x);
    int sy = (int)floor(y);

    // check bounds
    if (sx < 0 || sx >= width) return;
    if (sy < 0 || sy >= height) return;

    fill_pixel(sx, sy, color);
    return;
  }

  // Rasterize a line.
  void RasterizerImp::rasterize_line(float x0, float y0,
    float x1, float y1,
    Color color) {
    if (x0 > x1) {
      swap(x0, x1); swap(y0, y1);
    }

    float pt[] = { x0,y0 };
    float m = (y1 - y0) / (x1 - x0);
    float dpt[] = { 1,m };
    int steep = abs(m) > 1;
    if (steep) {
      dpt[0] = x1 == x0 ? 0 : 1 / abs(m);
      dpt[1] = x1 == x0 ? (y1 - y0) / abs(y1 - y0) : m / abs(m);
    }

    while (floor(pt[0]) <= floor(x1) && abs(pt[1] - y0) <= abs(y1 - y0)) {
      rasterize_point(pt[0], pt[1], color);
      pt[0] += dpt[0]; pt[1] += dpt[1];
    }
  }
  
  bool aboveLine(float x, float y, float x0, float y0, float x1, float y1) {
    return -(x - x0) * (y1 - y0) + (y - y0) * (x1 - x0) >= 0;
  }

  // Rasterize a triangle.
  void RasterizerImp::rasterize_triangle(float x0, float y0,
    float x1, float y1,
    float x2, float y2,
    Color color) {
    // TODO: Task 1: Implement basic triangle rasterization here, no supersampling
    
    if (!aboveLine(x2, y2, x0, y0, x1, y1)) {
      swap(x0, x2);
      swap(y0, y2);
    }
    
    int minX = (int)min(min(floor(x0), floor(x1)), floor(x2));
    int maxX = (int)max(max(ceil(x0), ceil(x1)), ceil(x2));
      
    int minY = (int)min(min(floor(y0), floor(y1)), floor(y2));
    int maxY = (int)max(max(ceil(y0), ceil(y1)), ceil(y2));
    
    int sampleGrid = (int) sqrt(sample_rate);
    float sampleDist = (float) 1 / sampleGrid;
    
    for (int x = max(minX, 0); x <= min(maxX, (int) width); x++) {
      for (int y = max(minY, 0); y <= min(maxY, (int) height); y++) {
//        std::cout << "test" << std::endl;
        
        for (int i = 0; i < sampleGrid; i++) {
          for (int j = 0; j < sampleGrid; j++) {
            float sampleX = x + sampleDist / 2 + sampleDist * i;
            float sampleY = y + sampleDist / 2 + sampleDist * j;
            if (!aboveLine(sampleX, sampleY, x0, y0, x1, y1)) {
              continue;
            }
            if (!aboveLine(sampleX, sampleY, x1, y1, x2, y2)) {
              continue;
            }
            if (!aboveLine(sampleX, sampleY, x2, y2, x0, y0)) {
              continue;
            }
        
            sample_buffer[sample_rate * (y * width + x) + i * sampleGrid + j] = color;
          }
        }
      }
    }
  }


  void RasterizerImp::rasterize_interpolated_color_triangle(float x0, float y0, Color c0,
    float x1, float y1, Color c1,
    float x2, float y2, Color c2)
  {
    // TODO: Task 4: Rasterize the triangle, calculating barycentric coordinates and using them to interpolate vertex colors across the triangle
    // Hint: You can reuse code from rasterize_triangle
    if (!aboveLine(x2, y2, x0, y0, x1, y1)) {
      swap(x0, x2);
      swap(y0, y2);
    }
    
    int minX = (int)min(min(floor(x0), floor(x1)), floor(x2));
    int maxX = (int)max(max(ceil(x0), ceil(x1)), ceil(x2));
      
    int minY = (int)min(min(floor(y0), floor(y1)), floor(y2));
    int maxY = (int)max(max(ceil(y0), ceil(y1)), ceil(y2));
    
    int sampleGrid = (int) sqrt(sample_rate);
    float sampleDist = (float) 1 / sampleGrid;
    
    for (int x = max(minX, 0); x <= min(maxX, (int) width); x++) {
      for (int y = max(minY, 0); y <= min(maxY, (int) height); y++) {
        for (int i = 0; i < sampleGrid; i++) {
          for (int j = 0; j < sampleGrid; j++) {
            float sampleX = x + sampleDist / 2 + sampleDist * i;
            float sampleY = y + sampleDist / 2 + sampleDist * j;
            if (!aboveLine(sampleX, sampleY, x0, y0, x1, y1)) {
              continue;
            }
            if (!aboveLine(sampleX, sampleY, x1, y1, x2, y2)) {
              continue;
            }
            if (!aboveLine(sampleX, sampleY, x2, y2, x0, y0)) {
              continue;
            }
            
            float a = (-(sampleX - x1) * (y2 - y1) + (sampleY - y1) * (x2 - x1)) / (-(x0 - x1) * (y2 - y1) + (y0 - y1) * (x2 - x1));
            
            float b = (-(sampleX - x2) * (y0 - y2) + (sampleY - y2) * (x0 - x2)) / (-(x1 - x2) * (y0 - y2) + (y1 - y2) * (x0 - x2));
            
            float c = 1 - a - b;

            sample_buffer[sample_rate * (y * width + x) + i * sampleGrid + j] = a * c0 + b * c1 + c * c2;
          }
        }
      }
    }
  }
  
  Vector3D findBarycentric(float x, float y, float x0, float y0, float x1, float y1, float x2, float y2) {
    float a = (-(x - x1) * (y2 - y1) + (y - y1) * (x2 - x1)) / (-(x0 - x1) * (y2 - y1) + (y0 - y1) * (x2 - x1));
    
    float b = (-(x - x2) * (y0 - y2) + (y - y2) * (x0 - x2)) / (-(x1 - x2) * (y0 - y2) + (y1 - y2) * (x0 - x2));
    
    float c = 1 - a - b;
    
    return Vector3D(a, b, c);
  }

  void RasterizerImp::rasterize_textured_triangle(float x0, float y0, float u0, float v0,
    float x1, float y1, float u1, float v1,
    float x2, float y2, float u2, float v2,
    Texture& tex)
  {
    // TODO: Task 5: Fill in the SampleParams struct and pass it to the tex.sample function.
    // TODO: Task 6: Set the correct barycentric differentials in the SampleParams struct.
    // Hint: You can reuse code from rasterize_triangle/rasterize_interpolated_color_triangle
    
    SampleParams sp;
    sp.psm = psm;
    sp.lsm = lsm;
    
    if (!aboveLine(x2, y2, x0, y0, x1, y1)) {
      swap(x0, x2);
      swap(y0, y2);
      swap(u0, u2);
      swap(v0, v2);
    }

    int minX = (int)min(min(floor(x0), floor(x1)), floor(x2));
    int maxX = (int)max(max(ceil(x0), ceil(x1)), ceil(x2));
      
    int minY = (int)min(min(floor(y0), floor(y1)), floor(y2));
    int maxY = (int)max(max(ceil(y0), ceil(y1)), ceil(y2));
    
    int sampleGrid = (int) sqrt(sample_rate);
    float sampleDist = (float) 1 / sampleGrid;
    
    for (int x = max(minX, 0); x < min(maxX, (int) width); x++) {
      for (int y = max(minY, 0); y < min(maxY, (int) height); y++) {
        for (int i = 0; i < sampleGrid; i++) {
          for (int j = 0; j < sampleGrid; j++) {
            float sampleX = x + sampleDist / 2 + sampleDist * i;
            float sampleY = y + sampleDist / 2 + sampleDist * j;
            if (!aboveLine(sampleX, sampleY, x0, y0, x1, y1)) {
              continue;
            }
            if (!aboveLine(sampleX, sampleY, x1, y1, x2, y2)) {
              continue;
            }
            if (!aboveLine(sampleX, sampleY, x2, y2, x0, y0)) {
              continue;
            }

            Vector3D b = findBarycentric(sampleX, sampleY, x0, y0, x1, y1, x2, y2);
            
            sp.p_uv = Vector2D(b[0] * u0 + b[1] * u1 + b[2] * u2, b[0] * v0 + b[1] * v1 + b[2] * v2);
            
            b = findBarycentric(sampleX + 1, sampleY, x0, y0, x1, y1, x2, y2);
            
            if (b[0] > 0 && b[1] > 0 && b[2] > 0) {
              sp.p_dx_uv = Vector2D(b[0] * u0 + b[1] * u1 + b[2] * u2, b[0] * v0 + b[1] * v1 + b[2] * v2) - sp.p_uv;
            }
            
            b = findBarycentric(sampleX, sampleY + 1, x0, y0, x1, y1, x2, y2);
            
            if (b[0] > 0 && b[1] > 0 && b[2] > 0) {
              sp.p_dy_uv = Vector2D(b[0] * u0 + b[1] * u1 + b[2] * u2, b[0] * v0 + b[1] * v1 + b[2] * v2) - sp.p_uv;
            }
            
            sample_buffer[sample_rate * (y * width + x) + i * sampleGrid + j] = tex.sample(sp);
          }
        }
      }
    }
  }

  void RasterizerImp::set_sample_rate(unsigned int rate) {
    // TODO: Task 2: You may want to update this function for supersampling support

    this->sample_rate = rate;
    clear_buffers();
    this->sample_buffer.resize(width * height * sample_rate, Color::White);
  }


  void RasterizerImp::set_framebuffer_target(unsigned char* rgb_framebuffer,
    size_t width, size_t height)
  {
    // TODO: Task 2: You may want to update this function for supersampling support

    this->width = width;
    this->height = height;
    this->rgb_framebuffer_target = rgb_framebuffer;
    
    clear_buffers();
    this->sample_buffer.resize(width * height * sample_rate, Color::White);
  }


  void RasterizerImp::clear_buffers() {
    std::fill(rgb_framebuffer_target, rgb_framebuffer_target + 3 * width * height, 255);
    std::fill(sample_buffer.begin(), sample_buffer.end(), Color::White);
  }


  // This function is called at the end of rasterizing all elements of the
  // SVG file.  If you use a supersample buffer to rasterize SVG elements
  // for antialising, you could use this call to fill the target framebuffer
  // pixels from the supersample buffer data.
  //
  void RasterizerImp::resolve_to_framebuffer() {
    // TODO: Task 2: You will likely want to update this function for supersampling support


    for (int x = 0; x < width; ++x) {
      for (int y = 0; y < height; ++y) {
        float r = 0;
        float g = 0;
        float b = 0;
    
        for (int i = 0; i < sample_rate; i++) {
          Color col = sample_buffer[sample_rate * (y * width + x) + i];
          r += col.r / sample_rate;
          g += col.g / sample_rate;
          b += col.b / sample_rate;
    
        }
    
        this->rgb_framebuffer_target[3 * (y * width + x)] = r * 255;
        this->rgb_framebuffer_target[3 * (y * width + x) + 1] = g * 255;
        this->rgb_framebuffer_target[3 * (y * width + x) + 2] = b * 255;
      }
    }
    


  }

  Rasterizer::~Rasterizer() { }


}// CGL
