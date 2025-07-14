import matplotlib.pyplot as plt
import numpy as np

x_coords = []
y_coords = []
with open('16_pts.txt', 'r') as f:
    for line in f:
        x, y = map(int, line.split())
        x_coords.append(x)
        y_coords.append(y)
x_coords = np.array(x_coords)
y_coords = np.array(y_coords)

# Create a scatter plot
plt.figure(figsize=(8, 8))

plt.scatter(x_coords, y_coords, c='blue', alpha=0.7, edgecolors='black')
# plt.scatter(x_coords * 16, y_coords * 16, c='red', alpha=0.7, edgecolors='black')
'''plt.scatter(x_coords * 256, y_coords * 256, c='green', alpha=0.7, edgecolors='black')
plt.scatter(x_coords * 4096, y_coords * 4096, c='yellow', alpha=0.7, edgecolors='black')
'''
# Add labels and title
plt.xlabel('X-axis')
plt.ylabel('Y-axis')

# Show grid for better visualization
plt.grid(True, linestyle='--', alpha=0.5)

# Display the plot
plt.show()
