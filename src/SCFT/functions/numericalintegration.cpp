#include "../SCFT.h"

double SCFT::trapz(double* f, int lower_bound, int upper_bound) {
  // If the two bounds are the same then return 0
  if(lower_bound == upper_bound) return 0;

  // Initialize
  int k_low = lower_bound;
  int k_up  = upper_bound;
  double out_val;

  // Add the start and end values
  out_val = 0.5*f[k_low] + 0.5*f[k_up];

  // Loop through the rest of the interval
  for(int k = k_low + 1; k < k_up; k++) {
      out_val += f[k];
  }

  // Return output
  return out_val;
}

double SCFT::trapz(double* z, double* f, int lower_bound, int upper_bound) {
  // If the two bounds are the same then return 0
  if(lower_bound == upper_bound) return 0;

  // Initialize
  int k_low = lower_bound;
  int k_up  = upper_bound;
  double out_val;

  // Initialize
  out_val = 0.0;

  // Loop through interval
  for(int k = k_low + 1; k <= k_up; k++) {
      out_val += (f[k] + f[k-1])*(z[k] - z[k-1])/2.0;
  }

  // Return output
  return out_val;
}

double SCFT::gaussQuad(double* z, double* f, int lower_bound, int upper_bound) {
  // If the two bounds are the same then return 0
  if(lower_bound == upper_bound) return 0;

  // Initialize
  int k_low = lower_bound;
  int k_up  = upper_bound;
  double out_val;

  // Initialize
  out_val = 0.0;

  // Loop through interval
  double low_xi = -1/sqrt(3);
  double high_xi = 1/sqrt(3);
  double temp;
  double df_low, df_high;
  double a,b,c,dz1,dz2;

  for(int k = k_low; k < k_up-1; k++) {
      // Calculate evaluation points
      low_xi = -(z[k+1] - z[k])/2.0/sqrt(3) + (z[k+1] + z[k])/2.0;
      high_xi = (z[k+1] - z[k])/2.0/sqrt(3) + (z[k+1] + z[k])/2.0;

      // Calculate slope at lower point and upper point
      if(k == 0) {
          dz1 = z[k+1] - z[k];
          dz2 = z[k+2] - z[k+1];
          a = -(2*dz1 + dz2)/dz1/(dz1 + dz2);
          b = (dz1 + dz2)/dz1/dz2;
          c = -dz1/dz2/(dz1 + dz2);
          df_low = a*f[k] + b*f[k+1] + c*f[k+2];

          dz1 = z[k+1] - z[k];
          dz2 = z[k+2] - z[k+1];
          a = -dz2/dz1/(dz1 + dz2);
          b = (dz2 - dz1)/dz1/dz2;
          c = dz1/dz2/(dz1 + dz2);
          df_high = a*f[k] + b*f[k+1] + c*f[k+2];
      }
      else if(k == k_up-2) {
          dz1 = z[k] - z[k-1];
          dz2 = z[k+1] - z[k];
          a = -dz2/dz1/(dz1 + dz2);
          b = (dz2 - dz1)/dz1/dz2;
          c = dz1/dz2/(dz1 + dz2);
          df_low = a*f[k-1] + b*f[k] + c*f[k+1];

          dz1 = z[k+1] - z[k];
          dz2 = z[k] - z[k-1];
          a = dz1/dz2/(dz1 + dz2);
          b = -(dz1 + dz2)/dz1/dz2;
          c = (2*dz1 + dz2)/dz1/(dz1 + dz2);
          df_high = a*f[k-1] + b*f[k] + c*f[k+1];
      }
      else {
          dz1 = z[k] - z[k-1];
          dz2 = z[k+1] - z[k];
          a = -dz2/dz1/(dz1 + dz2);
          b = (dz2 - dz1)/dz1/dz2;
          c = dz1/dz2/(dz1 + dz2);
          df_low = a*f[k-1] + b*f[k] + c*f[k+1];


          dz1 = z[k+1] - z[k];
          dz2 = z[k+2] - z[k+1];
          a = -dz2/dz1/(dz1 + dz2);
          b = (dz2 - dz1)/dz1/dz2;
          c = dz1/dz2/(dz1 + dz2);
          df_high = a*f[k] + b*f[k+1] + c*f[k+2];
      }

      // Interpolate to the two point quadrature
      temp = f[k] + df_low*(low_xi - z[k]);
      temp += f[k+1] + df_high*(high_xi - z[k+1]);

      // Add integral of sub-region to total
      out_val += (z[k+1] - z[k])/2.0*temp;
  }

  // Return output
  return out_val;
}


double SCFT::simps1(double* f, int lower_bound, int upper_bound) {
  // Adjust for small interval
  if(lower_bound == upper_bound) return 0;
  if(lower_bound - upper_bound == 1) return (f[lower_bound] + f[upper_bound])/2;

  // Declare variables
  int k_low, k_up, k_up0, oe_lower;
  bool odd_intervals;

  // Initialize the output
  double out_val;

  // Check the intervals
  oe_lower = lower_bound%2;
  odd_intervals = upper_bound%2 != lower_bound%2;

  // Limits of integration (Simpson's)
  k_up = upper_bound;
  k_up0 = k_up;
  k_low = lower_bound;

  // If odd number of intervals, then adjust bound
  if(odd_intervals) k_up0 = k_up - 1;

  // Add the start and end values
  out_val = f[k_low] + f[k_up0];

  // Loop through the rest of the interval
  for(int k = k_low + 1; k < k_up0; ++k) {
      out_val += f[k]*(2.0 + 2.0*(k%2 != oe_lower));
  }

  // Divide by 3 per Simpsons rule
  out_val /= 3.0;

  // If the final point was cut off to satisfy even intervals, trapezoid the last point
  if(odd_intervals)
  {
      out_val += (f[k_up] + f[k_up0])*0.5;
  }

  // Return output
  return out_val;
}

double SCFT::simps38(double* f, int lower_bound, int upper_bound) {
  // Adjust for small interval
  if(lower_bound == upper_bound) return 0;
  if(lower_bound - upper_bound == 1) return (f[lower_bound] + f[upper_bound])/2;
  if(lower_bound - upper_bound == 2) return (f[lower_bound] + 2.0*f[lower_bound+1] +f[upper_bound])/2;

  // Declare variables
  int k_low, k_up, k_up0;
  int divisThree;

  // Initialize the output
  double out_val;

  // Check the number of intervals
  divisThree = (upper_bound - lower_bound)%3;

  // Limits of integration
  k_up = upper_bound;
  k_up0 = k_up;
  k_low = lower_bound;

  // If number of intervals is not divisible by 3, then adjust bound
  k_up0 = k_up - divisThree;

  // Add the start and end values
  out_val = f[k_low] + f[k_up0];

  // Loop through the rest of the interval
  for(int k = k_low + 1; k < k_up0; k++) {
      out_val += f[k]*(3.0 - 1.0*((k-1)%3 == 2));
  }

  // Divide by 3 per Simpsons rule
  out_val *= 3.0/8.0;

  // We need to account for the last points if they were cut off. Use trapezoidal rule for this
  if(divisThree == 1) {
      out_val += (f[k_up0] + f[k_up])*0.5;
  }
  else if(divisThree == 2) {
      out_val += (f[k_up0] + 2*f[k_up0+1] + f[k_up])*0.5;
  }

  // Return output
  return out_val;
}
