//var macros = module.exports = {};

var macros = {
	upper: function () {
	  // `this` is the value in the parenthesis, or undefined if left out
	  return this.toUpperCase();
	},

	scale: function (percentage) {
	  var url = this;
	  return '<img src="' + url + '" style="width: ' + percentage + '" />';
	},

	color: function (color) {
	  return '<span style="color:' + color + '">' + this + '</span>'
	}
};
