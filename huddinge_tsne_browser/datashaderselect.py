import param
import holoviews as hv
from holoviews.streams import Stream
import logging as log


class DataShaderSelect(Stream):
    left_x, bottom_y, right_x, top_y = param.Number(), param.Number(
    ), param.Number(), param.Number()
    index = param.List(
        default=[],
        constant=True,
        doc="""
        Indices into a 1D datastructure.""")

    def __init__(self, source, dataset, x="tsne0", y="tsne1", **kwargs):
        super(DataShaderSelect, self).__init__(**kwargs)
        self._x = x
        self._y = y
        self._dataset = dataset
        self.__source = source
        self._source = hv.streams.Tap(source=source)
        self._source.add_subscriber(self.set_selection)
        self.selected = None
        log.info("Initialized DataShaderSelect")

    def set_selection(self, x, y):
        if x is not None and y is not None:
            d = self.__source.dframe()
            w_x = d.x.drop_duplicates().diff().max()
            w_y = d.y.drop_duplicates().diff().max()

            left_x, bottom_y, right_x, top_y = x - w_x / 2, y - w_y / 2, x + w_x / 2, y + w_y / 2
            q = "(x<={right_x})&(x>({left_x}))&(y<={top_y})&(y>({bottom_y}))".format(
                **locals())
            middle = d.query(q).iloc[0]
            self.left_x, self.bottom_y, self.right_x, self.top_y = middle.x - w_x / 2, middle.y - w_y / 2, middle.x + w_x / 2, middle.y + w_y / 2

            self.idx = (self._dataset[self._x] <= self.right_x) & (
                self._dataset[self._x] >= self.left_x) & (
                    self._dataset[self._y] <= self.top_y) & (
                        self._dataset[self._y] >= self.bottom_y)

        self.selected = list(self.idx.nonzero()[0])
        self.event(
            left_x=self.left_x,
            bottom_y=self.bottom_y,
            right_x=self.right_x,
            top_y=self.top_y,
            index=self.selected)
