class FancyListExtractions:
    """
    A utility class that provides extractions and manipulations for lists, including rolling slice assignments
    and finding indices of a subset within a main list.
    """

    @staticmethod
    def assign_rolling_slice(lst: list, start: int, new_values: list) -> list:
        """
        Assigns values from `new_values` to the list `lst` starting from the `start` index and rolls over
        if the end of the list is reached.

        Args:
                lst (list): The original list to modify.
                start (int): The starting index where the new values will be inserted.
                new_values (list): The list of new values to assign.

        Returns:
                list: The modified list with the new values assigned.

        Example::
        
                lst = [1, 2, 3, 4]
                new_values = [9, 10]
                result = FancyListExtractions.assign_rolling_slice(lst, start=3, new_values=new_values)
                # Output: [9, 2, 3, 10]
        """
        n = len(lst)
        length = len(new_values)

        # Loop through the new values and assign them to the original list in a rolling manner
        for i in range(length):
            lst[(start + i) % n] = new_values[i]

        return lst

    @staticmethod
    def find_indices(main_list: list, subset: list) -> list:
        """
        Finds the indices of the `subset` elements in the `main_list`. If a value in the subset is not found,
        returns `None` for that value's index.

        Args:
                main_list (list): The list to search within.
                subset (list): The list of values whose indices are to be found in the main list.

        Returns:
                list: A list of indices corresponding to the subset elements in the main list. If an element is not found, `None` is returned in its place.

        Example::
        
                main_list = ['a', 'b', 'c', 'd']
                subset = ['b', 'd', 'x']
                result = FancyListExtractions.find_indices(main_list, subset)
                # Output: [1, 3, None]
        """
        indices = []
        for value in subset:
            try:
                index = main_list.index(
                    value
                )  # Find the index of the value in the main list
                indices.append(index)
            except ValueError:
                indices.append(
                    None
                )  # If the value is not found, append None or handle as needed
        return indices
