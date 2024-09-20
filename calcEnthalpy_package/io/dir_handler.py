import os.path


class DirHandler:
	
	@staticmethod
	def mkdir_recursrive(folders, folder_path):
		
		if '/' != folder_path[-1]:
			folder_path += '/'
		if not os.path.exists(folder_path):
			os.mkdir(folder_path)
			
		for folder in folders:
			folder_path = folder_path + folder + '/'
			if not os.path.exists(folder_path):
				os.mkdir(folder_path)
			
		return folder_path
		