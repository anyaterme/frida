var createAllErrors = function() 
{
	var form = $( this ), errorList = $( "ul.errorMessages", form );

	var showAllErrorMessages = function() 
	{
		errorList.empty();
		// Find all invalid fields within the form.
		var invalidFields = form.find( ":invalid" ).each(function( index, node ) 
		{
			// Find the field's corresponding label
			var label = $( "label[for=" + node.id + "] "), // Opera incorrectly does not fill the validationMessage property.
			message = node.validationMessage || 'Invalid value.';
			reference = node.getAttribute("errorMsg");
			if (reference == "null")
				reference = ""

			errorList.show().append( "<li class='list-group-item list-group-item-danger'><span>" + label.html() + node.getAttribute("errorMsg") +  "</span> " + message + "</li>" );
		});
	};
	// Support Safari
	form.on( "submit", function( event ) 
	{
		if ( this.checkValidity && !this.checkValidity() ) 
		{
			$( this ).find( ":invalid" ).first().focus();
			event.preventDefault();
		}
	});
	$( "input[type=submit], button:not([type=button])", form).on( "click", showAllErrorMessages);
	$( "input", form ).on( "keypress", function( event ) 
	{
		var type = $( this ).attr( "type" );
		if ( /date|email|month|number|search|tel|text|time|url|week/.test(type) && event.keyCode == 13 ) 
		{
			showAllErrorMessages();
		}
	});
};

