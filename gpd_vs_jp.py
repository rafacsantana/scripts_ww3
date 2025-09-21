from scipy.stats import genpareto
import numpy as np

# Threshold
#u = np.percentile(X, 95)
#excesses = X[X > u] - u

# Fit GPD
#params = genpareto.fit(excesses)
#shape, loc, scale = params


import pandas as pd

# Load the data
file_path = "hatea_waiarohia.csv"
df = pd.read_csv(file_path)

# Display the first few rows and column names
df.head(), df.columns.tolist()

# Drop rows where any of the necessary columns (ss, ht, rain) are missing
df_clean = df.dropna(subset=["ss", "ht", "rain"]).copy()

# Compute sea level as skew surge + high tide
df_clean["sea_level"] = df_clean["ss"] + df_clean["ht"]

# Keep only the relevant variables
df_extremes = df_clean[["year", "month", "day", "sea_level", "rain"]].copy()

# Show summary and sample
df_extremes.describe(), df_extremes.head()

#######################################################################

#import numpy as np
#from scipy.stats import genpareto

# Define sea level variable
sea_level = df_extremes["sea_level"].values

# Use 95th percentile as threshold
threshold = np.percentile(sea_level, 95)

# Extract exceedances above the threshold
exceedances = sea_level[sea_level > threshold] - threshold

# Fit GPD to exceedances
gpd_shape, gpd_loc, gpd_scale = genpareto.fit(exceedances)

# Store threshold and parameters
gpd_params = {
    "threshold": threshold,
    "shape": gpd_shape,
    "scale": gpd_scale,
    "loc": gpd_loc
}

print(gpd_params)


########################################################################

# Calculate the 0.1% AEP (Annual Exceedance Probability) return level from GPD
# AEP = 0.001 → Return period = 1 / 0.001 = 1000-year
# N = number of data points
N = len(sea_level)
exceedance_count = len(exceedances)
prob_exceed = exceedance_count / N  # Empirical exceedance probability above threshold

# Return level formula for GPD:
# z_p = u + (σ/ξ) * ((p / exceed_prob)^(-ξ) - 1)  if ξ ≠ 0
p = 0.001  # 0.1% AEP

# Use the formula
if gpd_shape != 0:
    z_p = threshold + (gpd_scale / gpd_shape) * ((p / prob_exceed) ** (-gpd_shape) - 1)
else:
    z_p = threshold - gpd_scale * np.log(p / prob_exceed)

print(z_p)



import matplotlib.pyplot as plt
from sklearn.linear_model import LinearRegression
#import scipy.stats.linregress as LinearRegression
import seaborn as sns

# Filter data for conditional model: sea_level > threshold
df_extreme = df_clean[df_clean["sea_level"] > threshold].copy()

# Log-transform rainfall + small constant to stabilize variance
df_extreme["log_rain"] = np.log1p(df_extreme["rain"])

# Fit linear model: log(rain) ~ sea_level (approximation to conditional extremes model)
X = df_extreme["sea_level"].values.reshape(-1, 1)
y = df_extreme["log_rain"].values
reg = LinearRegression().fit(X, y)

# Predict expected rainfall at 0.1% AEP sea level
predicted_log_rain = reg.predict([[return_level_gpd]])
predicted_rain = np.expm1(predicted_log_rain[0])  # inverse of log1p

# Prepare joint return curve data (for plotting)
# Simulate sea levels across a range and predict associated rainfall
sea_range = np.linspace(threshold, return_level_gpd + 0.1, 100).reshape(-1, 1)
rain_log_preds = reg.predict(sea_range)
rain_preds = np.expm1(rain_log_preds)

# Save joint curve data for plotting
joint_curve_df = pd.DataFrame({
    "sea_level": sea_range.flatten(),
    "expected_rain": rain_preds
})

# Plot joint return curve
plt.figure(figsize=(8, 6))
sns.scatterplot(data=df_extreme, x="sea_level", y="rain", alpha=0.3, label="Observed")
plt.plot(joint_curve_df["sea_level"], joint_curve_df["expected_rain"],
         color="red", lw=2, label="Joint Return Curve (approx.)")
plt.scatter([return_level_gpd], [predicted_rain], color="black", label="0.1% AEP Point")
plt.xlabel("Sea Level (m)")
plt.ylabel("Rainfall (mm)")
plt.title("Joint Return Curve: Rainfall vs Sea Level")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

predicted_rain


