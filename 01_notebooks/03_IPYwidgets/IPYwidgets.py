#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# #Manual and interactive segmentation
# import ipywidgets as widgets
# from IPython.display import display
# import numpy as np
# import matplotlib.pyplot as plt
# import cv2
# import tifffile

class ManualSegmentation():
 
    def __init__(self, image, cmap='Spectral', polygon_color=(255, 0, 0)):
        self.ipython = get_ipython()
        if self.ipython:
            self.ipython.run_line_magic('matplotlib', 'widget')

        self.image = self.process_image(image)
        self.polygon_color = polygon_color
        self.selected_points = []

        self.figure_to_draw_points, self.axes_in_figure = plt.subplots(figsize=(5, 5))
        self.new_image = self.axes_in_figure.imshow(self.image, cmap=cmap)

        self.click = self.figure_to_draw_points.canvas.mpl_connect('button_press_event', self.onclick)

    def process_image(self, image):
#         processed_image = RemoveExtrema(image, min_percentile=0.1, max_percentile=99).remove_outliers()
        processed_image = image
        processed_image = (processed_image - processed_image.min()) / (processed_image.max() - processed_image.min()) * 255
        return processed_image

    def polygon(self, new_image, points_in_polygon):
        points_in_polygon = np.array(points_in_polygon, np.int32)
        points_in_polygon = points_in_polygon.reshape((-1, 1, 2))
        cv2.polylines(new_image, [points_in_polygon], isClosed=True, color=self.polygon_color, thickness=3)
        return new_image

    def switch_to_inline(self):
        if self.ipython:
            self.ipython.run_line_magic('matplotlib', 'inline')
        plt.show()  # Ensure that any existing plots are displayed
        return None

    def onclick(self, event):
        if event.xdata is not None and event.ydata is not None:
            self.selected_points.append([int(event.xdata), int(event.ydata)])
            updated_image = self.polygon(np.copy(self.image), self.selected_points)
            
            for point in self.selected_points:
                cv2.circle(updated_image, tuple(point), radius=3, color=self.polygon_color, thickness=-1)
            
            self.new_image.set_data(updated_image)
        return None

    def close_and_save(self, filename='temp_mask.tif', save_mask=True):
        if self.selected_points:
            # Create an empty array with the same shape as the image slice
            mask = np.zeros(self.image.shape[:2], dtype=np.uint8)

            # Convert selected points to a numpy array for cv2.fillPoly
            mask_array = np.array([self.selected_points], dtype=np.int32)

            # Create the mask
            cv2.fillPoly(mask, mask_array, 255)  # Fill with white (255)

            if save_mask:
                # Save the mask as a tif file
                tifffile.imwrite(filename, mask, dtype='uint8')
                print(f'Mask saved as {filename}')

            # Close the figure and disconnect the click
            self.figure_to_draw_points.canvas.mpl_disconnect(self.click)
            plt.close(self.figure_to_draw_points)
            plt.close()

            self.switch_to_inline()  # Switch back to inline when closing the plot
            return mask.astype(bool)
        else:
            print('No points selected to create a mask.')
            plt.close()
            return None


# Example usage:
mask_object = ManualSegmentation(
    image=np.max(image_colors[2,:,:,:],axis=0), 
    cmap='Spectral'
)  # NOTICE THAT FOR THE MASK SELECTION WE USE THE GREEN CHANNEL (channel 1)


# # #Test
# import ipywidgets as widgets
# from IPython.display import display

# # Create a simple slider widget
# slider = widgets.IntSlider(value=50, min=0, max=100, step=1, description='Slider:')

# # Display the slider
# display(slider)

# # Function to update based on slider value
# def on_value_change(change):
#     print(f'Slider value changed to: {change["new"]}')

# # Attach the function to the slider
# slider.observe(on_value_change, names='value')

